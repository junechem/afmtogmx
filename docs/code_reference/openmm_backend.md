# `openmm_backend.py` — `OpenMMBackend` (`off.openmm`)

The public OpenMM API. Accessed as `off.openmm` on a `ReadOFF`. Owns OpenMM config;
reads shared FF data from `self._parent.{bonded, nonbonded, charges}`. Delegates to
`xml_generation.py` (XML build) and `pdb_processing.py` (PDB prep).

## Construction / state
`__init__(self, parent)` — stores the parent and initializes `self.config`:
- `incl_mol=[]` — molecule names to include (empty = all).
- `molname_translations={}` — maps `.off` molecule names to PDB residue names,
  e.g. `{'H2OQM': 'SOL'}`.

Same explicit→config→default parameter resolution pattern as `GROMACSBackend`.

## Config methods
- `set_config(**kwargs)` — updates config in place; returns `self`.
- `get_config(key=None)` — one value, or a copy of the whole config dict.

## Methods
- `gen_xml(output='forcefield.xml', incl_mol=None, molname_translations=None)` —
  builds a complete OpenMM `<ForceField>` XML from parsed `.off` data and writes it.
  Atom types are namespaced `"<MOLNAME>_<TYPE>"` to avoid collisions. Section order:
  AtomTypes → Residues → NonbondedForce → Bond/Angle/Dihedral forces → EXP/SRD/STR
  custom nonbonded forces. Calls the `xml_generation.*` collectors and builders, then
  `write_xml`.
  - Supported nonbonded types: EXP, STR/STRC, SRD, POW (folded into SRD with r0=0),
    BUC (split into EXP + SRD).
  - Units converted: kcal/mol→kJ/mol (×4.184), Å→nm (×0.1).
  - Charges live on `parent.charges` keyed by atom **name**; mapped name→type via
    the ATO section before writing `<NonbondedForce>`.
- `preprocess_pdb(input_pdb, xml_file, output_pdb='output.pdb', maxwarn=0)` —
  renames PDB atoms to match the XML residue atom names and emits fresh CONECT
  records. Reads **only** the XML (no in-memory state), so one XML can process many
  PDBs. Returns `self` (chainable). Just two calls:
  `pdb_processing.build_residue_topology_from_xml` then `pdb_processing.process_pdb`.
  - Atom identification strategy is chosen automatically: **topology matching**
    (graph isomorphism) when the PDB has CONECT records — also reorders atoms into
    canonical XML order; **positional matching** (Nth atom → Nth XML name) otherwise,
    with a `maxwarn` element-mismatch failsafe.
  - `maxwarn` only applies on the positional path. Raises `FileNotFoundError`,
    `IOError`, or `SystemExit` (unmatchable residue / too many positional mismatches).
