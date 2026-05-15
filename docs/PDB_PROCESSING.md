# PDB Processing

This document explains everything the `afmtogmx` package can do with PDB
files. Read this instead of re-deriving the behavior from the source.

## Purpose

OpenMM simulations that use custom forces need a PDB file that:

1. Uses **atom names that match the force-field XML's `<Residues>` section**, and
2. Contains **explicit `CONECT` records** describing the covalent topology.

Standard PDB files coming out of packing/equilibration tools usually have
neither. PDB processing rewrites such a file so it is OpenMM-ready, using an
OpenMM-style force-field XML as the single source of truth.

## Where the code lives

| File | Role |
|------|------|
| `src/afmtogmx/core/pdb_processing.py` | The actual implementation — stateless module-level functions, no classes. |
| `src/afmtogmx/core/openmm_backend.py` | Thin wrapper: `OpenMMBackend.preprocess_pdb()`, reached as `off.openmm.preprocess_pdb(...)`. |

> Note: older standalone copies of a class-based `pdb_preprocessor.py` exist
> under `test/openmm*/FF_input/`. Those are superseded artifacts. The
> supported implementation is `pdb_processing.py`.

## The two-stage pipeline

PDB processing is deliberately split into two independent stages so a single
parsed XML can be reused across many PDB files.

### Stage 1 — `build_residue_topology_from_xml(xml_file)`

Parses an OpenMM-style `<ForceField>` XML and returns a per-residue topology
dict:

```python
{resname: {'atom_names': [ordered_names], 'bonds': [(n1, n2), ...]}}
```

- `atom_names` — every `<Atom name=...>` inside `<Residues>/<Residue>`, **in XML order**.
- `bonds` — every `<Bond atomName1=... atomName2=.../>` child of that `<Residue>`.

### Stage 2 — `process_pdb(input_pdb, output_pdb, residue_topology)`

Consumes the topology dict from Stage 1 and rewrites the PDB:

1. **Renames atoms positionally.** For each residue instance, the *Nth* atom
   in the PDB receives the *Nth* name from `atom_names`. There is **no
   name-based matching** — order is everything.
2. **Strips** all existing `TER`, `CONECT`, and `END` records.
3. **Regenerates `CONECT` records** from the topology's `bonds`.
4. **Writes the output** in this order: non-atom header records (TITLE,
   REMARK, CRYST1, MODEL, ...) → ATOM/HETATM records → CONECT records → `END`.

## Public API

### Via a `ReadOFF` object (recommended)

```python
import afmtogmx as afm

off = afm.ReadOFF(off_loc='butanol_water.off')
off.load_charges_from_file('charges.txt')
off.openmm.set_config(molname_translations={'H2OQM': 'SOL'})
off.openmm.gen_xml(output='forcefield.xml')

# Rename atoms + add CONECT records using the XML just generated.
off.openmm.preprocess_pdb('conf.pdb', 'forcefield.xml', 'conf_processed.pdb')
```

`preprocess_pdb(input_pdb, xml_file, output_pdb='output.pdb')`:
- `input_pdb` — PDB to process.
- `xml_file` — an OpenMM force-field XML (typically produced by `off.openmm.gen_xml()`, but any OpenMM-style XML works).
- `output_pdb` — destination path (default `'output.pdb'`).
- Returns `self` for method chaining.
- Internally just calls Stage 1 then Stage 2 — it consults **no** in-memory `ReadOFF` state.

### Standalone (no `ReadOFF` needed)

Because the functions read everything from the XML, they can be used on their
own whenever a force-field XML is available:

```python
from afmtogmx.core import pdb_processing

topology = pdb_processing.build_residue_topology_from_xml('forcefield.xml')
pdb_processing.process_pdb('conf.pdb', 'conf_processed.pdb', topology)
```

This is the intended path when generating many PDBs against one XML — parse
the XML once, reuse `topology`.

## Key design decisions (and why)

**Connectivity comes from `<Residue>`'s own `<Bond>` children — not from the
force sections.** An OpenMM XML has bond information in two places:

- `<Residues>/<Residue>/<Bond atomName1 atomName2>` — actual covalent topology, by atom *name*.
- `<HarmonicBondForce>` / `<CustomBondForce>` `<Bond>` entries — force *parameters*, by atom *type*/*class*.

PDB processing uses the **residue-level atom-name bonds**. Using the
type-pair force entries as connectivity would cause **Cartesian-product bond
explosion** for residues where several atoms share one type (e.g. all H atoms
typed the same).

**Positional renaming, not name matching.** Stage 2 assumes the PDB lists a
residue's atoms in the same order as the XML `<Residue>` block. This keeps the
logic simple and is reliable for AFM-generated systems where both files come
from the same molecule definition.

**Residue instances are keyed by `(residue_name, chain_id, residue_seq)`.**
Each instance gets its own position counter and its own `CONECT` set, so
repeated copies of the same residue are handled independently.

**`CONECT` records are deduplicated globally** by a `frozenset` of the two
atom indices, so a bond defined twice yields a single `CONECT` line. The final
list is sorted before writing.

## PDB column parsing

`process_pdb` reads fixed-width PDB columns from `ATOM`/`HETATM` lines:

| Field | Columns (0-indexed slice) |
|-------|---------------------------|
| Atom serial number | `[6:11]` |
| Atom name | `[12:16]` |
| Residue name | `[17:20]` |
| Chain ID | `[21:22]` |
| Residue sequence number | `[22:26]` |

When renaming, the new name is written right-justified into a 4-wide field at
columns `[12:16]` (`f"{name:>4}"`); all other columns are preserved verbatim.

## Behavior on edge cases / limitations

- **Residue not present in the XML topology** — its atoms are left with their
  original names and produce **no** `CONECT` records. The residue is passed
  through untouched. (Note: unlike the old standalone script, the current
  module prints no warning for this.)
- **More PDB atoms in a residue than the XML defines** — once the position
  counter passes the end of `atom_names`, the extra atoms keep their original
  names. No error is raised.
- **PDB atom order differs from XML residue order** — atoms are silently
  renamed to the *wrong* names. Stage 2 cannot detect this; the input must be
  correctly ordered.
- **No atom class→type conversion.** The current module does positional
  renaming only. (An earlier standalone preprocessor attempted class-to-type
  conversion; that logic was dropped.)
- **Virtual sites (mass = 0)** are treated like any other atom — they are
  counted in the positional renaming.
- Existing `TER`/`CONECT`/`END` records are always discarded and regenerated.

## Internal functions reference

All in `pdb_processing.py`:

| Function | Public? | Role |
|----------|---------|------|
| `build_residue_topology_from_xml(xml_file)` | yes | Stage 1 orchestrator — XML → topology dict. |
| `process_pdb(input_pdb, output_pdb, residue_topology)` | yes | Stage 2 — rewrite the PDB. |
| `_parse_residues_from_xml(root)` | no | Parse `<Residues>` into `{resname: {atom_names, bonds}}`. |
| `_parse_bond_parameters_from_xml(root)` | no | Parse `<HarmonicBondForce>`/`<CustomBondForce>` type-pair bonds. **Currently unused** — kept as a seam for a future XML→`ReadOFF` parser that also rebuilds parameters. |
| `_generate_conect_records(atom_records, residue_topology)` | no | Build the sorted, deduplicated `CONECT` lines. |

The module is intentionally one-helper-per-XML-section so that a future
expansion (parsing angles, dihedrals, nonbonded forces, etc.) can add helpers
alongside the existing ones without disturbing the orchestrator.

## End-to-end workflow context

PDB processing is the second half of the `.off` → OpenMM workflow:

1. `off.openmm.gen_xml(...)` — convert parsed `.off` force-field data into an OpenMM `<ForceField>` XML.
2. `off.openmm.preprocess_pdb(...)` — rename atoms + add `CONECT` records in the simulation PDB so it matches that XML.

After both steps, the XML and the processed PDB can be loaded directly into an
OpenMM simulation that uses the custom AFM force field.
