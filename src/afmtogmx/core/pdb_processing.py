"""PDB preprocessing for OpenMM compatibility.

Module-level helpers that prepare a PDB file for OpenMM simulation by:
1. Renaming atoms to match force-field XML atom definitions.
2. Adding CONECT records from force-field bond topology.

This module reads everything from an OpenMM-style XML force field. It has
no dependency on :class:`ReadOFF` or any other in-memory force-field
object, so it can be used standalone whenever a force-field XML file is
available — typical when generating many PDBs against a single XML.

XML parsing is split into one helper per section
(:func:`_parse_residues_from_xml`, :func:`_parse_bonds_from_xml`, ...) so
that a future expansion to read the entire force field (angles,
dihedrals, nonbonded, etc.) can drop in alongside the existing helpers
without disturbing the orchestrator.
"""

import xml.etree.ElementTree as ET


# --------------------------------------------------------------------------
# XML section parsers — one per top-level XML section
# --------------------------------------------------------------------------

def _parse_residues_from_xml(root):
    """Parse the ``<Residues>`` section.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        Root element of the parsed force-field XML.

    Returns
    -------
    dict
        ``{resname: {'atom_names': [...], 'type_to_names': {type: [names]}}}``
        ``atom_names`` preserves XML order. ``type_to_names`` allows mapping
        a bond's type-pair back to actual atom names within the residue.
    """
    residues = {}
    for residue in root.findall('Residues/Residue'):
        resname = residue.attrib['name']
        atom_names = []
        type_to_names = {}
        for atom in residue.findall('Atom'):
            aname = atom.attrib['name']
            atype = atom.attrib['type']
            atom_names.append(aname)
            type_to_names.setdefault(atype, []).append(aname)
        residues[resname] = {
            'atom_names': atom_names,
            'type_to_names': type_to_names,
        }
    return residues


def _parse_bonds_from_xml(root):
    """Parse all bond-defining force sections into a list of type pairs.

    Handles both attribute conventions:
    - ``<HarmonicBondForce><Bond class1=... class2=.../></HarmonicBondForce>``
    - ``<CustomBondForce><Bond type1=... type2=.../></CustomBondForce>``

    The ``<Bond>`` child elements inside ``<Residue>`` are *not* read here;
    those reference atom names rather than types and are stored separately
    by OpenMM's force-field definition.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element

    Returns
    -------
    list of tuple
        ``[(type1, type2), ...]`` — each entry is a pair of force-field
        atom type names that participate in a bond.
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

    Joins the residue definitions (atom names + types) with the global
    bond list (type-pair → bond) so that each residue ends up with the
    explicit list of (atom_name, atom_name) bonds within it.

    Parameters
    ----------
    xml_file : str
        Path to an OpenMM-style ``<ForceField>`` XML file.

    Returns
    -------
    dict
        ``{resname: {'atom_names': [ordered_names], 'bonds': [(n1, n2), ...]}}``
        in the format expected by :func:`process_pdb`.

    Examples
    --------
    >>> topology = build_residue_topology_from_xml('forcefield.xml')
    >>> process_pdb('conf.pdb', 'conf_processed.pdb', topology)
    """
    root = ET.parse(xml_file).getroot()
    residues = _parse_residues_from_xml(root)
    bond_type_pairs = _parse_bonds_from_xml(root)

    result = {}
    for resname, info in residues.items():
        type_to_names = info['type_to_names']
        bonds = []
        for t1, t2 in bond_type_pairs:
            for n1 in type_to_names.get(t1, []):
                for n2 in type_to_names.get(t2, []):
                    if n1 != n2:
                        bonds.append((n1, n2))
        result[resname] = {'atom_names': info['atom_names'], 'bonds': bonds}
    return result


# --------------------------------------------------------------------------
# PDB rewriter — consumes a residue_topology dict, emits the new PDB
# --------------------------------------------------------------------------

def process_pdb(input_pdb, output_pdb, residue_topology):
    """Rename PDB atoms to match residue topology and generate CONECT records.

    Atoms are renamed *positionally*: the Nth atom of each residue instance
    in the input PDB receives the Nth name from
    ``residue_topology[resname]['atom_names']``. Existing TER/CONECT/END
    records are stripped; new CONECT records are emitted from the topology.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file.
    output_pdb : str
        Path for output PDB file.
    residue_topology : dict
        From :func:`build_residue_topology_from_xml`.

    Raises
    ------
    FileNotFoundError
        If ``input_pdb`` does not exist.
    IOError
        If ``output_pdb`` cannot be written.
    """
    atom_records = []
    other_records = []
    residue_atom_positions = {}

    with open(input_pdb, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            atom_index = int(line[6:11].strip())
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21:22]
            residue_seq = line[22:26].strip()

            converted_atom_name = atom_name
            residue_key = (residue_name, chain_id, residue_seq)
            position = residue_atom_positions.setdefault(residue_key, 0)

            if residue_name in residue_topology:
                xml_atom_list = residue_topology[residue_name]['atom_names']
                if position < len(xml_atom_list):
                    converted_atom_name = xml_atom_list[position]
                    line = line[:12] + f"{converted_atom_name:>4}" + line[16:]
                residue_atom_positions[residue_key] += 1

            atom_records.append({
                'line': line,
                'index': atom_index,
                'name': converted_atom_name,
                'residue': residue_name,
                'chain': chain_id,
                'seq': residue_seq,
            })
        elif not line.startswith(('TER', 'CONECT', 'END')):
            other_records.append(line)

    conect_records = _generate_conect_records(atom_records, residue_topology)

    with open(output_pdb, 'w') as f:
        for record in other_records:
            f.write(record)
        for atom in atom_records:
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
