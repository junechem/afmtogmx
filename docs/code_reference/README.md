# Code Reference — `src/afmtogmx/core/`

These files document every function in the `core/` modules so a future reader can
understand the codebase **without re-reading every source file**. Each reference
file mirrors one `.py` file and lists, per function: its signature, what it does,
its inputs/outputs, and any non-obvious behavior.

> **Maintenance rule:** When you change a function in `core/`, update its entry in
> the matching reference file in the same commit. When you add or remove a
> function, add or remove its entry. These docs are only useful if they stay in
> sync with the code. If a doc and the code disagree, the code is the source of
> truth — fix the doc.

## Module map

| Source file | Reference | Role |
|-------------|-----------|------|
| `gen_md.py` | [gen_md.md](gen_md.md) | Main `ReadOFF` class: parse `.off`, hold FF data, dispatch to backends |
| `functions.py` | [functions.md](functions.md) | Low-level `.off` parsing + interaction math (exp/srd/shtr/powe/quarbond) |
| `tabulated_potentials.py` | [tabulated_potentials.md](tabulated_potentials.md) | GROMACS tabulated potential (`.xvg`) generation, bonded + nonbonded |
| `topology.py` | [topology.md](topology.md) | GROMACS `.top` topology file generation |
| `gmx_backend.py` | [gmx_backend.md](gmx_backend.md) | `GROMACSBackend` (`off.gmx`): config + public GROMACS API |
| `openmm_backend.py` | [openmm_backend.md](openmm_backend.md) | `OpenMMBackend` (`off.openmm`): config + public OpenMM API |
| `xml_generation.py` | [xml_generation.md](xml_generation.md) | Build OpenMM `<ForceField>` XML sections from FF data |
| `pdb_processing.py` | [pdb_processing.md](pdb_processing.md) | Preprocess PDB for OpenMM: rename atoms, emit CONECT |
| `populate_bonded.py` | [populate_bonded.md](populate_bonded.md) | Build a new molecule's bonded dict from a user topology directory |
| `report.py` | [report.md](report.md) | Publication-style fixed-width text export of the parsed force field |

Not documented here (excluded by request / trivial): `compare.py`, `residues.py`,
`__init__.py` (empty).

## How the pieces fit together

```
ReadOFF(off_loc)                      # gen_md.py — parses .off via functions.py
  ├── off.bonded / off.nonbonded      # parsed FF data (shared by both backends)
  ├── off.charges / off.residues
  ├── off.gmx   -> GROMACSBackend     # gmx_backend.py
  │      ├── gen_nonbonded_tabpot / gen_bonded_tabpot   -> tabulated_potentials.py
  │      ├── write_*_tabpot                             -> tabulated_potentials.py
  │      └── gen_nonbonded_topology / gen_bonded_topology -> topology.py
  ├── off.openmm -> OpenMMBackend     # openmm_backend.py
  │      ├── gen_xml          -> xml_generation.py
  │      └── preprocess_pdb   -> pdb_processing.py
  └── off.write_report                -> report.py
```

`ReadOFF.change_molecule` and `ReadOFF.populate_bonded` (the latter via
`populate_bonded.py`) rewrite the shared FF data before either backend runs.
`ReadOFF.write_report` is engine-independent — it renders the parsed data as text
rather than producing simulation input.
