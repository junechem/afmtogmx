# `pdb_processing.py` — PDB preprocessing for OpenMM

Module-level helpers (no class) that prepare a PDB so OpenMM can apply a
generated force-field XML: they **rename PDB atoms** to the XML residue atom
names and **emit fresh CONECT records** from the XML bond topology. Driven by
`OpenMMBackend.preprocess_pdb` (`off.openmm`), which calls
`build_residue_topology_from_xml` then `process_pdb`.

**Standalone:** depends only on `xml.etree.ElementTree` — no `ReadOFF`, no
OpenMM import. Reads everything it needs from the XML file, so one XML can
process many PDBs.

## Two atom-identification strategies (chosen automatically in `process_pdb`)
- **Topology matching** — used when the PDB has CONECT records. Each residue's
  bond graph is matched to the XML residue graph by element-labelled graph
  isomorphism, identifying atoms by *connectivity* (robust to differing atom
  order, e.g. CIF-derived PDBs); atoms are then reordered into canonical XML
  order.
- **Positional matching** — fallback when there are no CONECT records. Nth PDB
  atom → Nth XML name, guarded by a GROMACS-style `maxwarn` element-mismatch
  failsafe.

## Element helpers
- `_ELEMENT_SYMBOLS` — module-level `frozenset` of all real element symbols.
  Kept local so the module needs no periodic-table dependency.
- `_norm_element(symbol)` — `symbol.strip().capitalize()` (e.g. `'CL'` → `'Cl'`).
- `_element_from_name(atom_name)` — element from an atom **name**. Takes the
  leading alphabetic characters, then resolves to a real element by trying the
  first **two** characters as a symbol before falling back to the first
  character. Keeps two-letter elements intact (`'NA'`→`'Na'`, `'CL'`→`'Cl'`)
  while still mapping multi-letter names whose prefix is a single-letter element
  (`'OW'`→`'O'`, `'HW'`→`'H'`, `'O0'`→`'O'`). A name with no valid element prefix
  (e.g. the virtual site `'EW'`) falls back to the capitalized first character
  (`'E'`); empty name → `''`. **Used identically on the PDB side and the XML
  virtual-site fallback so both agree.**
- `_atom_element(line)` — element of a PDB ATOM/HETATM line: reads the element
  column (77–78); if blank, falls back to `_element_from_name` on the atom name.

## XML section parsers (one per top-level XML section)
- `_parse_atomtypes_from_xml(root)` — `<AtomTypes>` → `{type_name: element}`.
  Virtual-site types (no `element` attribute) are omitted.
- `_parse_residues_from_xml(root, type_to_element)` — `<Residues>` →
  `{resname: {'atom_names': [...], 'bonds': [(n1, n2), ...], 'elements':
  {name: element}}}`. Connectivity comes from each `<Residue>`'s own
  `<Bond atomName1=… atomName2=…/>` children (the actual covalent topology) —
  **not** the type-pair `<Bond>` entries in the force sections, which carry only
  parameters and would cause Cartesian-product bond explosion. Per-atom element
  resolves via `type_to_element`, falling back to `_element_from_name`.
- `_parse_bond_parameters_from_xml(root)` — `<HarmonicBondForce>` /
  `<CustomBondForce>` `<Bond>` type-pairs → `[(type1, type2), ...]`. **Not** used
  for connectivity; a seam for a future XML→`ReadOFF` parameter parser.

## Orchestrator
- `build_residue_topology_from_xml(xml_file)` — parse the XML and return the
  `residue_topology` dict (shape above) for `process_pdb`. Just
  `_parse_atomtypes_from_xml` then `_parse_residues_from_xml`.

## PDB record parsing
- `_parse_atom_line(line)` — ATOM/HETATM line → record dict
  (`line, index, name, residue, chain, seq, element`).
- `_parse_conect_records(lines)` — all CONECT lines → undirected edge set of
  `frozenset({serial_i, serial_j})` (first serial is the central atom).

## Graph isomorphism
- `_match_residue_graph(pdb_atoms, pdb_edges, xml_entry)` — element-labelled
  graph isomorphism by recursive backtracking (VF2-lite). XML atoms tried in
  descending-degree order, paired with an unused PDB atom of equal element and
  degree whose adjacency is consistent with all prior pairings. Returns
  `{pdb_serial: xml_name}` (first complete mapping) or `None`. Equivalent atoms
  share a type, so any valid isomorphism is acceptable.

## Atom-name assignment strategies
- `_rename(atom, new_name)` — rewrite the name field (columns 13–16) in place.
- `_assign_names_by_topology(residues, residue_order, conect_edges,
  residue_topology)` — CONECT path. Renames in place via `_match_residue_graph`.
  `exit(1)` if a residue's atom count disagrees with the XML or no isomorphism
  is found (with a diagnostic naming the residue and bond counts).
- `_assign_names_by_position(residues, residue_order, residue_topology,
  maxwarn)` — no-CONECT path. Nth atom → Nth XML name; per residue **type**,
  counts atoms whose element disagrees with the assigned name. Each mismatched
  residue type is one warning; if warning count exceeds `maxwarn`, prints a
  detailed failure and `exit(1)`, otherwise prints warnings and continues.

## PDB rewriter (top-level)
- `process_pdb(input_pdb, output_pdb, residue_topology, maxwarn=0)` — the main
  routine. Parses ATOM/HETATM and CONECT, groups atoms by residue instance
  (first-seen order), then:
  - Warns for each PDB residue **absent** from the XML (atoms passed through
    unprocessed); `exit(1)` if **no** PDB residue matches any XML residue
    (usually a wrong name / missing `molname_translations`).
  - Picks topology vs positional path by presence of CONECT edges.
  - Reorders topology-matched residues into canonical XML order; other residues
    keep original order. Renumbers serials sequentially.
  - Strips existing TER/CONECT/END, writes non-atom records, atoms, freshly
    generated CONECT records, then `END`.
  - Raises `FileNotFoundError` / `IOError` / `SystemExit` (see above).
- `_generate_conect_records(atom_records, residue_topology)` — emit CONECT lines
  for every bond of every residue instance, deduplicated globally by
  `frozenset((idx1, idx2))`; returns them sorted.
