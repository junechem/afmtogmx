"""PDB preprocessing for OpenMM compatibility.

Module-level helpers that prepare a PDB file for OpenMM simulation by:
1. Renaming atoms to match force-field XML atom definitions.
2. Adding CONECT records from force-field bond topology.

This module reads everything from an OpenMM-style XML force field. It has
no dependency on :class:`ReadOFF` or any other in-memory force-field
object, so it can be used standalone whenever a force-field XML file is
available — typical when generating many PDBs against a single XML.

:func:`process_pdb` chooses between two atom-identification strategies:

* **Topology matching** — used when the input PDB carries CONECT
  records. Each residue's bond graph is matched against the XML residue
  graph by element-labelled graph isomorphism, so atoms are identified
  by *connectivity* rather than by their order in the file. Atoms are
  then reordered into the canonical XML order. This is robust to PDBs
  whose atom order differs from the ``.off``/XML order (e.g. PDBs built
  from a CIF).
* **Positional matching** — fallback used when the PDB has no CONECT
  records. The Nth atom of each residue is given the Nth XML name. A
  GROMACS-style ``maxwarn`` failsafe compares each atom's element with
  the element of the name it would receive and aborts on mismatch
  unless enough warnings are permitted.

XML parsing is split into one helper per section so that a future
expansion to read the entire force field (angles, dihedrals, nonbonded,
etc.) can drop in alongside the existing helpers.
"""

import xml.etree.ElementTree as ET


# --------------------------------------------------------------------------
# Element helpers
# --------------------------------------------------------------------------

def _norm_element(symbol):
    """Normalise an element symbol for comparison (``'CL'`` → ``'Cl'``)."""
    return symbol.strip().capitalize()


def _element_from_name(atom_name):
    """Best-effort element guess from an atom name (``'O0'`` → ``'O'``)."""
    element = ''
    for ch in atom_name:
        if ch.isalpha():
            element += ch
        else:
            break
    return _norm_element(element)


def _atom_element(line):
    """Element of a PDB ATOM/HETATM line.

    Reads the element column (77-78); falls back to the atom name when
    that column is blank or absent.
    """
    element = line[76:78].strip() if len(line) > 76 else ''
    if element:
        return _norm_element(element)
    return _element_from_name(line[12:16].strip())


# --------------------------------------------------------------------------
# XML section parsers — one per top-level XML section
# --------------------------------------------------------------------------

def _parse_atomtypes_from_xml(root):
    """Parse the ``<AtomTypes>`` section into ``{type_name: element}``.

    Virtual-site types carry no ``element`` attribute and are simply
    omitted from the result.
    """
    type_to_element = {}
    for atom_type in root.findall('AtomTypes/Type'):
        name = atom_type.attrib.get('name')
        element = atom_type.attrib.get('element')
        if name and element:
            type_to_element[name] = _norm_element(element)
    return type_to_element


def _parse_residues_from_xml(root, type_to_element):
    """Parse the ``<Residues>`` section.

    Reads each ``<Residue>``'s atom list (preserving XML order), its
    per-atom element (via the atom's ``type`` and the ``<AtomTypes>``
    table, falling back to the name), and its per-residue
    ``<Bond atomName1=... atomName2=.../>`` connectivity. These
    atom-name bonds are the actual covalent topology — distinct from the
    type-pair entries inside ``<HarmonicBondForce>`` / ``<CustomBondForce>``,
    which carry only the bonded-force parameters and would generate
    Cartesian-product bond explosion if used as connectivity for residues
    with multiple atoms sharing a type.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        Root element of the parsed force-field XML.
    type_to_element : dict
        ``{type_name: element}`` from :func:`_parse_atomtypes_from_xml`.

    Returns
    -------
    dict
        ``{resname: {'atom_names': [...], 'bonds': [(n1, n2), ...],
        'elements': {name: element}}}``
    """
    residues = {}
    for residue in root.findall('Residues/Residue'):
        resname = residue.attrib['name']
        atom_names = []
        elements = {}
        for atom in residue.findall('Atom'):
            aname = atom.attrib['name']
            atom_names.append(aname)
            atype = atom.attrib.get('type')
            elements[aname] = (type_to_element.get(atype)
                               or _element_from_name(aname))
        bonds = [
            (bond.attrib['atomName1'], bond.attrib['atomName2'])
            for bond in residue.findall('Bond')
        ]
        residues[resname] = {
            'atom_names': atom_names,
            'bonds': bonds,
            'elements': elements,
        }
    return residues


def _parse_bond_parameters_from_xml(root):
    """Parse the global bond *parameter* sections into a list of type pairs.

    Reads ``<HarmonicBondForce>``/``<CustomBondForce>`` ``<Bond>`` entries,
    which assign force-field parameters to atom-type pairs. NOT used for
    PDB connectivity (see :func:`_parse_residues_from_xml`); kept as a
    seam for a future XML→ReadOFF parser that also rebuilds parameters.

    Returns
    -------
    list of tuple
        ``[(type1, type2), ...]``
    """
    bonds = []
    for bond in root.findall('HarmonicBondForce/Bond'):
        bonds.append((bond.attrib['class1'], bond.attrib['class2']))
    for bond in root.findall('CustomBondForce/Bond'):
        bonds.append((bond.attrib['type1'], bond.attrib['type2']))
    return bonds


# --------------------------------------------------------------------------
# Orchestrator — combines parsed sections into the residue_topology shape
# --------------------------------------------------------------------------

def build_residue_topology_from_xml(xml_file):
    """Read an OpenMM XML force field and return per-residue topology.

    The result is suitable as input to :func:`process_pdb`. Connectivity
    is taken straight from each ``<Residue>``'s own ``<Bond>`` children
    (atom-name bonds), which is the canonical place OpenMM stores residue
    topology — independent of any force-section parameter assignments.
    Per-atom elements are resolved through the ``<AtomTypes>`` table.

    Parameters
    ----------
    xml_file : str
        Path to an OpenMM-style ``<ForceField>`` XML file.

    Returns
    -------
    dict
        ``{resname: {'atom_names': [...], 'bonds': [(n1, n2), ...],
        'elements': {name: element}}}``

    Examples
    --------
    >>> topology = build_residue_topology_from_xml('forcefield.xml')
    >>> process_pdb('conf.pdb', 'conf_processed.pdb', topology)
    """
    root = ET.parse(xml_file).getroot()
    type_to_element = _parse_atomtypes_from_xml(root)
    return _parse_residues_from_xml(root, type_to_element)


# --------------------------------------------------------------------------
# PDB record parsing
# --------------------------------------------------------------------------

def _parse_atom_line(line):
    """Parse a PDB ATOM/HETATM line into a record dict."""
    return {
        'line': line,
        'index': int(line[6:11].strip()),
        'name': line[12:16].strip(),
        'residue': line[17:20].strip(),
        'chain': line[21:22],
        'seq': line[22:26].strip(),
        'element': _atom_element(line),
    }


def _parse_conect_records(lines):
    """Parse all CONECT lines into an undirected edge set.

    Returns a set of ``frozenset({serial_i, serial_j})``. The first
    serial on a CONECT line is the central atom; every remaining serial
    is one of its bonded neighbours.
    """
    edges = set()
    for line in lines:
        if not line.startswith('CONECT'):
            continue
        try:
            serials = [int(tok) for tok in line[6:].split()]
        except ValueError:
            continue
        for neighbour in serials[1:]:
            if neighbour != serials[0]:
                edges.add(frozenset((serials[0], neighbour)))
    return edges


# --------------------------------------------------------------------------
# Graph isomorphism — identify atoms by bond topology
# --------------------------------------------------------------------------

def _match_residue_graph(pdb_atoms, pdb_edges, xml_entry):
    """Match a residue's PDB bond graph to its XML graph by isomorphism.

    Element-labelled graph isomorphism via recursive backtracking
    (VF2-lite). XML atoms are tried in descending-degree order; each is
    paired with an unused PDB atom of equal element and degree whose
    adjacency is consistent with every pair already matched. The first
    complete mapping found is returned — equivalent atoms (e.g. the
    hydrogens on one carbon) share a force-field type, so any valid
    isomorphism is acceptable.

    Parameters
    ----------
    pdb_atoms : list of dict
        Atom records (need ``'index'`` and ``'element'``) for one residue.
    pdb_edges : set of frozenset
        Bonds among those atoms, as ``frozenset({serial_i, serial_j})``.
    xml_entry : dict
        A ``residue_topology`` entry (``atom_names``, ``bonds``, ``elements``).

    Returns
    -------
    dict or None
        ``{pdb_serial: xml_name}`` mapping, or ``None`` if no isomorphism
        exists.
    """
    xml_names = xml_entry['atom_names']
    xml_elements = xml_entry['elements']
    xml_adj = {n: set() for n in xml_names}
    for n1, n2 in xml_entry['bonds']:
        if n1 in xml_adj and n2 in xml_adj:
            xml_adj[n1].add(n2)
            xml_adj[n2].add(n1)

    pdb_indices = [a['index'] for a in pdb_atoms]
    pdb_element = {a['index']: a['element'] for a in pdb_atoms}
    pdb_adj = {i: set() for i in pdb_indices}
    for edge in pdb_edges:
        i, j = tuple(edge)
        if i in pdb_adj and j in pdb_adj:
            pdb_adj[i].add(j)
            pdb_adj[j].add(i)

    order = sorted(xml_names, key=lambda n: -len(xml_adj[n]))
    mapping = {}   # xml_name -> pdb_index
    used = set()   # pdb indices already assigned

    def backtrack(pos):
        if pos == len(order):
            return True
        name = order[pos]
        for idx in pdb_indices:
            if idx in used:
                continue
            if pdb_element[idx] != _norm_element(xml_elements[name]):
                continue
            if len(pdb_adj[idx]) != len(xml_adj[name]):
                continue
            consistent = True
            for mapped_name, mapped_idx in mapping.items():
                xml_adjacent = mapped_name in xml_adj[name]
                pdb_adjacent = mapped_idx in pdb_adj[idx]
                if xml_adjacent != pdb_adjacent:
                    consistent = False
                    break
            if consistent:
                mapping[name] = idx
                used.add(idx)
                if backtrack(pos + 1):
                    return True
                del mapping[name]
                used.discard(idx)
        return False

    if backtrack(0):
        return {idx: name for name, idx in mapping.items()}
    return None


# --------------------------------------------------------------------------
# Atom-name assignment strategies
# --------------------------------------------------------------------------

def _rename(atom, new_name):
    """Rewrite an atom record's name field (columns 13-16) in place."""
    atom['line'] = atom['line'][:12] + f"{new_name:>4}" + atom['line'][16:]
    atom['name'] = new_name


def _assign_names_by_topology(residues, residue_order, conect_edges,
                              residue_topology):
    """Identify atoms by bond-graph isomorphism (CONECT-based path).

    Renames atoms in place. Aborts with a clear message if a residue's
    atom count disagrees with the XML or no isomorphism is found.
    """
    for key in residue_order:
        resname, chain, seq = key
        if resname not in residue_topology:
            continue
        atoms = residues[key]
        xml_entry = residue_topology[resname]
        n_xml = len(xml_entry['atom_names'])

        if len(atoms) != n_xml:
            print(f"PDB preprocessing failed: residue {resname} "
                  f"(chain {chain!r}, seq {seq}) has {len(atoms)} atom(s) "
                  f"but the force-field XML defines {n_xml}.")
            exit(1)

        indices = {a['index'] for a in atoms}
        pdb_edges = {e for e in conect_edges if e <= indices}
        mapping = _match_residue_graph(atoms, pdb_edges, xml_entry)

        if mapping is None:
            print(f"PDB preprocessing failed: could not match residue "
                  f"{resname} (chain {chain!r}, seq {seq}) to the force-field "
                  f"XML by bond topology. The residue has {len(pdb_edges)} "
                  f"CONECT bond(s); the XML defines {len(xml_entry['bonds'])}. "
                  f"Check that the PDB CONECT records describe the same "
                  f"molecule as the XML.")
            exit(1)

        for atom in atoms:
            _rename(atom, mapping[atom['index']])


def _assign_names_by_position(residues, residue_order, residue_topology,
                              maxwarn):
    """Identify atoms by position (no-CONECT fallback path).

    The Nth atom of each residue receives the Nth XML name. Each atom's
    element is compared with the element of its assigned name; a residue
    *type* with any mismatch produces one warning. If the warning count
    exceeds ``maxwarn`` the run aborts, otherwise warnings are printed
    and processing continues.
    """
    stats = {}  # resname -> {'mismatch': int, 'atoms': int, 'instances': int}

    for key in residue_order:
        resname = key[0]
        if resname not in residue_topology:
            continue
        atoms = residues[key]
        xml_names = residue_topology[resname]['atom_names']
        xml_elements = residue_topology[resname]['elements']

        s = stats.setdefault(resname,
                              {'mismatch': 0, 'atoms': 0, 'instances': 0})
        s['instances'] += 1

        for pos, atom in enumerate(atoms):
            if pos >= len(xml_names):
                break
            new_name = xml_names[pos]
            s['atoms'] += 1
            if atom['element'] != _norm_element(xml_elements[new_name]):
                s['mismatch'] += 1
            _rename(atom, new_name)

    warnings = [
        (f"Residue {resname}: {s['mismatch']} of {s['atoms']} positional "
         f"atom-name assignments (across {s['instances']} residue instance(s)) "
         f"have an element that disagrees with the XML")
        for resname, s in stats.items() if s['mismatch'] > 0
    ]

    if len(warnings) > maxwarn:
        message = [
            "PDB preprocessing failed: positional atom-name assignment looks "
            "wrong.",
            f"  {len(warnings)} warning(s), maxwarn={maxwarn}:",
        ]
        message += [f"    - {w}" for w in warnings]
        message += [
            "  The input PDB has no CONECT records, so atoms were matched by",
            "  position against the force-field XML; the element mismatches",
            "  above mean the PDB atom order likely differs from the XML.",
            "  Fix the atom ordering, add CONECT records to enable topology",
            f"  matching, or pass maxwarn>={len(warnings)} to override.",
        ]
        print('\n'.join(message))
        exit(1)

    for w in warnings:
        print(f"Warning: {w}")


# --------------------------------------------------------------------------
# PDB rewriter — orchestrates parsing, naming, reordering, CONECT output
# --------------------------------------------------------------------------

def process_pdb(input_pdb, output_pdb, residue_topology, maxwarn=0):
    """Rename PDB atoms to match residue topology and generate CONECT records.

    When the input PDB carries CONECT records, atoms are identified by
    bond-graph isomorphism against the XML residue graph and each residue
    is reordered into the canonical XML atom order. When it does not,
    atoms are renamed *positionally* (Nth atom → Nth XML name) and a
    ``maxwarn`` element-mismatch failsafe guards against mis-ordered
    input. Existing TER/CONECT/END records are stripped; new CONECT
    records are emitted from the topology and serial numbers are
    renumbered sequentially.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file.
    output_pdb : str
        Path for output PDB file.
    residue_topology : dict
        From :func:`build_residue_topology_from_xml`.
    maxwarn : int, optional
        Number of element-mismatch warnings tolerated on the positional
        (no-CONECT) path before the run aborts. Default ``0`` (strict).
        Ignored when the PDB has CONECT records.

    Raises
    ------
    FileNotFoundError
        If ``input_pdb`` does not exist.
    IOError
        If ``output_pdb`` cannot be written.
    SystemExit
        If a residue cannot be matched to the XML (topology path) or the
        positional element-mismatch count exceeds ``maxwarn``.
    """
    with open(input_pdb, 'r') as f:
        lines = f.readlines()

    atom_records = []
    other_records = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            atom_records.append(_parse_atom_line(line))
        elif not line.startswith(('TER', 'CONECT', 'END')):
            other_records.append(line)

    conect_edges = _parse_conect_records(lines)

    # Group atoms by residue instance, preserving first-seen order.
    residues = {}
    residue_order = []
    for atom in atom_records:
        key = (atom['residue'], atom['chain'], atom['seq'])
        if key not in residues:
            residues[key] = []
            residue_order.append(key)
        residues[key].append(atom)

    # Guard against a residue-name mismatch silently producing a no-op:
    # warn for each PDB residue absent from the XML, and abort outright if
    # nothing at all matches (almost always a wrong name / missing
    # molname_translations).
    topology_resnames = set(residue_topology)
    pdb_resnames = {key[0] for key in residue_order}
    matched_instances = sum(1 for key in residue_order
                            if key[0] in topology_resnames)

    for resname in sorted(pdb_resnames - topology_resnames):
        print(f"Warning: PDB residue {resname!r} is not defined in the "
              f"force-field XML; its atoms are passed through unprocessed "
              f"(XML residues: {', '.join(sorted(topology_resnames))}).")

    if matched_instances == 0:
        print("PDB preprocessing failed: no PDB residue matches any residue "
              "in the force-field XML.\n"
              f"  PDB residues: {', '.join(sorted(pdb_resnames))}\n"
              f"  XML residues: {', '.join(sorted(topology_resnames))}\n"
              "  Check the residue names, or set molname_translations so the "
              "XML residue name matches the PDB.")
        exit(1)

    use_topology = bool(conect_edges)
    if use_topology:
        _assign_names_by_topology(residues, residue_order, conect_edges,
                                  residue_topology)
    else:
        _assign_names_by_position(residues, residue_order, residue_topology,
                                  maxwarn)

    # Build the final atom list: topology-matched residues are reordered
    # into canonical XML order; everything else keeps its original order.
    final_atoms = []
    for key in residue_order:
        resname = key[0]
        atoms = residues[key]
        if use_topology and resname in residue_topology:
            by_name = {a['name']: a for a in atoms}
            final_atoms.extend(by_name[name]
                               for name in residue_topology[resname]['atom_names'])
        else:
            final_atoms.extend(atoms)

    # Renumber serials sequentially so CONECT records stay consistent.
    for new_index, atom in enumerate(final_atoms, start=1):
        atom['index'] = new_index
        atom['line'] = atom['line'][:6] + f"{new_index:>5}" + atom['line'][11:]

    conect_records = _generate_conect_records(final_atoms, residue_topology)

    with open(output_pdb, 'w') as f:
        for record in other_records:
            f.write(record)
        for atom in final_atoms:
            f.write(atom['line'])
        for conect in conect_records:
            f.write(conect)
        f.write("END\n")


def _generate_conect_records(atom_records, residue_topology):
    """Emit CONECT lines for every bond in every residue instance.

    Bonds are deduplicated globally by ``frozenset`` of (idx1, idx2) so a
    bond defined twice in the topology (e.g. via duplicate entries in two
    bond sections) yields one CONECT line.
    """
    conect_records = []
    bonds_added = set()

    residues = {}
    for atom in atom_records:
        res_key = (atom['residue'], atom['chain'], atom['seq'])
        residues.setdefault(res_key, []).append(atom)

    for (residue_name, _chain, _seq), atoms in residues.items():
        if residue_name not in residue_topology:
            continue
        name_to_index = {atom['name']: atom['index'] for atom in atoms}
        for n1, n2 in residue_topology[residue_name]['bonds']:
            if n1 in name_to_index and n2 in name_to_index:
                idx1, idx2 = name_to_index[n1], name_to_index[n2]
                bond_key = frozenset((idx1, idx2))
                if bond_key not in bonds_added:
                    bonds_added.add(bond_key)
                    conect_records.append(f"CONECT{idx1:5}{idx2:5}\n")

    return sorted(conect_records)
