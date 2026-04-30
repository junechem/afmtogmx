"""PDB preprocessing for OpenMM compatibility.

Module-level helpers that prepare a PDB file for OpenMM simulation by:
1. Renaming atoms to match force-field XML atom definitions.
2. Adding CONECT records from force-field bond topology.

Used by :class:`afmtogmx.core.openmm_backend.OpenMMBackend`.
"""


def _unique_atom_names(atom_map):
    """Per-residue unique names like ``C0, C1, H0, H1`` (element + counter).

    Matches the naming scheme used in xml_generation._unique_atom_names.
    """
    counters, result = {}, {}
    for atid in sorted(atom_map.keys()):
        atname, attype = atom_map[atid]
        if attype in ('NETF', 'TORQ'):
            continue
        element = attype[0]
        n = counters.get(element, 0)
        counters[element] = n + 1
        result[atid] = f'{element}{n}'
    return result


def build_residue_topology(bonded, mol_names, molname_translations):
    """Build atom name ordering and bond topology from bonded force-field data.

    Returns a dict mapping residue names to atom ordering and bonds, derived
    directly from bonded structures rather than by re-parsing XML. This ensures
    consistency with gen_xml() output.

    Parameters
    ----------
    bonded : dict
        Dictionary from ReadOFF.bonded with structure
        {mol_name: {ATO, BON, ...}}.
    mol_names : list of str
        Molecule names to include (filtered by incl_mol in caller).
    molname_translations : dict
        Maps .off molecule names to PDB-compatible residue names.

    Returns
    -------
    dict
        {resname: {'atom_names': [ordered_names], 'bonds': [(at1, at2), ...]}}
        Keys are translated residue names; atom_names follow PDB column layout
        (element + counter); bonds are tuples of unique atom names.
    """
    result = {}

    for mol in mol_names:
        resname = molname_translations.get(mol, mol)
        atom_map = {int(k): v for k, v in bonded[mol]['ATO']['All'].items()}
        unique = _unique_atom_names(atom_map)

        atom_names = [unique[atid] for atid in sorted(unique.keys())]

        bonds = []
        for bond_type, pairs_dict in bonded[mol].get('BON', {}).items():
            for params, atom_pairs in pairs_dict.items():
                for at1_id, at2_id in atom_pairs:
                    if at1_id not in unique or at2_id not in unique:
                        continue
                    bonds.append((unique[at1_id], unique[at2_id]))

        result[resname] = {'atom_names': atom_names, 'bonds': bonds}

    return result


def process_pdb(input_pdb, output_pdb, residue_topology):
    """Rename PDB atoms to match residue topology and generate CONECT records.

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file.
    output_pdb : str
        Path for output PDB file.
    residue_topology : dict
        Output from build_residue_topology().

    Returns
    -------
    None
        Writes output PDB file in place.
    """
    atom_records = []
    other_records = []
    residue_atom_positions = {}

    try:
        with open(input_pdb, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"PDB file not found: {input_pdb}")

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            record_type = line[0:6].strip()
            atom_index = int(line[6:11].strip())
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21:22]
            residue_seq = line[22:26].strip()

            converted_atom_name = atom_name

            residue_key = (residue_name, chain_id, residue_seq)
            if residue_key not in residue_atom_positions:
                residue_atom_positions[residue_key] = 0

            if residue_name in residue_topology:
                xml_atom_list = residue_topology[residue_name]['atom_names']
                position = residue_atom_positions[residue_key]

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
                'seq': residue_seq
            })
        elif not line.startswith(('TER', 'CONECT', 'END')):
            other_records.append(line)

    conect_records = _generate_conect_records(atom_records, residue_topology)

    try:
        with open(output_pdb, 'w') as f:
            for record in other_records:
                if not record.startswith('END'):
                    f.write(record)
            for atom in atom_records:
                f.write(atom['line'])
            for conect in conect_records:
                f.write(conect)
            f.write("END\n")
    except IOError as e:
        raise IOError(f"Error writing output file: {e}")


def _generate_conect_records(atom_records, residue_topology):
    """Generate CONECT records from atom records and residue bond definitions.

    Parameters
    ----------
    atom_records : list of dict
        Atom records with keys: index, name, residue, chain, seq.
    residue_topology : dict
        Output from build_residue_topology().

    Returns
    -------
    list of str
        Sorted CONECT record strings (one per bond).
    """
    conect_records = []
    bonds_added = set()

    residues = {}
    for atom in atom_records:
        res_key = (atom['residue'], atom['chain'], atom['seq'])
        if res_key not in residues:
            residues[res_key] = []
        residues[res_key].append(atom)

    for (residue_name, chain, seq), atoms in residues.items():
        if residue_name not in residue_topology:
            continue

        xml_atom_to_pdb_index = {}
        for atom in atoms:
            atom_name = atom['name']
            pdb_index = atom['index']
            xml_atom_to_pdb_index[atom_name] = pdb_index

        for atom1_name, atom2_name in residue_topology[residue_name]['bonds']:
            if atom1_name in xml_atom_to_pdb_index and atom2_name in xml_atom_to_pdb_index:
                idx1 = xml_atom_to_pdb_index[atom1_name]
                idx2 = xml_atom_to_pdb_index[atom2_name]

                bond_key = tuple(sorted([idx1, idx2]))
                if bond_key not in bonds_added:
                    conect_records.append(f"CONECT{idx1:5}{idx2:5}\n")
                    bonds_added.add(bond_key)

    return sorted(conect_records)
