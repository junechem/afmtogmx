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

    Reads each ``<Residue>``'s atom list (preserving XML order) and its
    per-residue ``<Bond atomName1=... atomName2=.../>`` connectivity.
    These atom-name bonds are the actual covalent topology — distinct
    from the type-pair entries inside ``<HarmonicBondForce>`` /
    ``<CustomBondForce>``, which carry only the bonded-force parameters
    and would generate Cartesian-product bond explosion if used as
    connectivity for residues with multiple atoms sharing a type.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        Root element of the parsed force-field XML.

    Returns
    -------
    dict
        ``{resname: {'atom_names': [...], 'bonds': [(n1, n2), ...]}}``
    """
    residues = {}
    for residue in root.findall('Residues/Residue'):
        resname = residue.attrib['name']
        atom_names = [atom.attrib['name'] for atom in residue.findall('Atom')]
        bonds = [
            (bond.attrib['atomName1'], bond.attrib['atomName2'])
            for bond in residue.findall('Bond')
        ]
        residues[resname] = {'atom_names': atom_names, 'bonds': bonds}
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

    Parameters
    ----------
    xml_file : str
        Path to an OpenMM-style ``<ForceField>`` XML file.

    Returns
    -------
    dict
        ``{resname: {'atom_names': [ordered_names], 'bonds': [(n1, n2), ...]}}``

    Examples
    --------
    >>> topology = build_residue_topology_from_xml('forcefield.xml')
    >>> process_pdb('conf.pdb', 'conf_processed.pdb', topology)
    """
    root = ET.parse(xml_file).getroot()
    return _parse_residues_from_xml(root)


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
