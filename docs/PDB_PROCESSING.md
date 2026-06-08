# PDB Processing

This document explains everything the `afmtogmx` package can do with PDB
files. Read this instead of re-deriving the behavior from the source.

## Purpose

OpenMM simulations that use custom forces need a PDB file that:

1. Uses **atom names that match the force-field XML's `<Residues>` section**, and
2. Contains **explicit `CONECT` records** describing the covalent topology.

Standard PDB files coming out of packing/equilibration tools — or built from a
CIF — usually have neither, and often list atoms in a different order than the
`.off`/XML. PDB processing rewrites such a file so it is OpenMM-ready, using an
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

PDB processing is split into two independent stages so a single parsed XML can
be reused across many PDB files.

### Stage 1 — `build_residue_topology_from_xml(xml_file)`

Parses an OpenMM-style `<ForceField>` XML and returns a per-residue topology
dict:

```python
{resname: {'atom_names': [ordered_names],
           'bonds': [(n1, n2), ...],
           'elements': {name: element}}}
```

- `atom_names` — every `<Atom name=...>` inside `<Residues>/<Residue>`, **in XML order**.
- `bonds` — every `<Bond atomName1=... atomName2=.../>` child of that `<Residue>`.
- `elements` — element of each atom, resolved through the atom's `type` and the
  `<AtomTypes>` table (`<Type element=.../>`), falling back to the leading
  letters of the atom name.

### Stage 2 — `process_pdb(input_pdb, output_pdb, residue_topology, maxwarn=0)`

Consumes the topology dict from Stage 1 and rewrites the PDB. It picks one of
**two atom-identification strategies automatically**, based on whether the
input PDB has `CONECT` records:

**Path A — topology matching (CONECT records present).**
For each residue instance, the bond graph is read from the `CONECT` records and
matched against the XML residue graph by **element-labelled graph
isomorphism**. This identifies each physical atom by its *connectivity*, so it
works no matter what order the atoms appear in. Atoms are then **reordered into
the canonical XML order**.

**Path B — positional matching (no CONECT records).**
The *Nth* atom of each residue receives the *Nth* name from `atom_names`. A
GROMACS-style `maxwarn` failsafe (see below) guards against mis-ordered input.
The original atom order is kept.

Both paths then:

1. **Strip** all existing `TER`, `CONECT`, and `END` records.
2. **Renumber** atom serial numbers sequentially (`1..N`).
3. **Regenerate `CONECT` records** from the topology's `bonds`.
4. **Write the output** in this order: non-atom header records (TITLE, REMARK,
   CRYST1, MODEL, ...) → ATOM/HETATM records → CONECT records → `END`.

### The `maxwarn` failsafe (Path B only)

On the positional path, each atom's element is compared with the element of the
XML name it would receive. A residue *type* with any mismatch produces **one
warning**. If the number of warnings exceeds `maxwarn` (default `0`),
`process_pdb` aborts with a message that lists the mismatches and explains the
likely cause (PDB atom order ≠ XML order). Raising `maxwarn` lets a
genuinely-correct-but-unverifiable PDB through — much like `grompp -maxwarn`.

`maxwarn` is **ignored on the topology path**: graph matching already verifies
elements as part of the isomorphism.

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

`preprocess_pdb(input_pdb, xml_file, output_pdb='output.pdb', maxwarn=0)`:
- `input_pdb` — PDB to process.
- `xml_file` — an OpenMM force-field XML (typically produced by `off.openmm.gen_xml()`, but any OpenMM-style XML works).
- `output_pdb` — destination path (default `'output.pdb'`).
- `maxwarn` — element-mismatch warnings tolerated on the positional path (default `0`; ignored when the PDB has CONECT records).
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

**Topology matching when CONECT exists; positional only as a fallback.** A PDB
built from a CIF can list atoms in any order — and may even use a *different*
order in each residue instance. Positional renaming would then silently assign
names to the wrong atoms. When `CONECT` records are present they give the exact
bond graph, which uniquely identifies every atom; matching that graph against
the XML residue graph is order-independent and correct. Positional matching is
kept only for PDBs that lack `CONECT`, and is guarded by `maxwarn`.

**Graph isomorphism is element-labelled and degree-pruned.** `_match_residue_graph`
is a small VF2-style backtracking matcher: XML atoms are tried in
descending-degree order and paired with PDB atoms of equal element and degree
whose adjacency is consistent with every pair already matched. Molecules are
tiny (tens of atoms), so this is effectively instant even for thousands of
residue instances. Chemically equivalent atoms (e.g. the hydrogens on one
carbon) yield several valid isomorphisms; the first one found is used, which is
safe because equivalent atoms share a force-field type.

**Residue instances are keyed by `(residue_name, chain_id, residue_seq)`.**
Each instance is matched, reordered, and given its own `CONECT` set
independently, so repeated copies of the same residue are handled separately.

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
| Element symbol | `[76:78]` (falls back to the atom name if blank) |

`CONECT` lines are parsed by whitespace-splitting the serials after the
`CONECT` tag; the first serial is the central atom, the rest are its bonded
neighbours.

When renaming, the new name is written right-justified into a 4-wide field at
columns `[12:16]`; the renumbered serial into `[6:11]`; all other columns are
preserved verbatim.

## Behavior on edge cases / limitations

- **Residue not present in the XML topology** — its atoms are left with their
  original names, are not reordered, and produce **no** `CONECT` records. The
  residue is passed through (serials are still renumbered). No warning is
  printed.
- **Topology path, atom-count mismatch** — if a residue instance has a
  different number of atoms than the XML residue, `process_pdb` aborts
  (`SystemExit`) naming the residue.
- **Topology path, no isomorphism** — if a residue's `CONECT` graph cannot be
  matched to the XML residue graph (e.g. wrong/missing bonds), `process_pdb`
  aborts naming the residue and reporting the bond counts on each side. This is
  a hard error — not `maxwarn`-gated.
- **Positional path, element mismatch** — counted per residue type; aborts when
  warnings exceed `maxwarn`. See "The `maxwarn` failsafe" above.
- **Positional path, more PDB atoms than the XML defines** — atoms past the end
  of `atom_names` keep their original names; no error is raised.
- **No atom class→type conversion.** The module assigns atom *names*; the XML
  residue maps name → type, and OpenMM resolves the rest.
- **Virtual sites (mass = 0)** are treated like any other atom — included in
  matching, renaming, and reordering.
- Existing `TER`/`CONECT`/`END` records are always discarded and regenerated.

## Internal functions reference

All in `pdb_processing.py`:

| Function | Public? | Role |
|----------|---------|------|
| `build_residue_topology_from_xml(xml_file)` | yes | Stage 1 orchestrator — XML → topology dict. |
| `process_pdb(input_pdb, output_pdb, residue_topology, maxwarn=0)` | yes | Stage 2 — rewrite the PDB. |
| `_parse_atomtypes_from_xml(root)` | no | Parse `<AtomTypes>` into `{type_name: element}`. |
| `_parse_residues_from_xml(root, type_to_element)` | no | Parse `<Residues>` into `{resname: {atom_names, bonds, elements}}`. |
| `_parse_bond_parameters_from_xml(root)` | no | Parse `<HarmonicBondForce>`/`<CustomBondForce>` type-pair bonds. **Currently unused** — a seam for a future XML→`ReadOFF` parser that also rebuilds parameters. |
| `_parse_atom_line(line)` | no | Parse one ATOM/HETATM line into a record dict. |
| `_parse_conect_records(lines)` | no | Parse all `CONECT` lines into an undirected edge set. |
| `_match_residue_graph(pdb_atoms, pdb_edges, xml_entry)` | no | Element-labelled graph isomorphism — `{pdb_serial: xml_name}` or `None`. |
| `_assign_names_by_topology(...)` | no | Path A — rename via graph matching; abort on mismatch. |
| `_assign_names_by_position(...)` | no | Path B — positional rename + `maxwarn` element failsafe. |
| `_generate_conect_records(atom_records, residue_topology)` | no | Build the sorted, deduplicated `CONECT` lines. |
| `_norm_element` / `_element_from_name` / `_atom_element` / `_rename` | no | Small element/record helpers. |

The module is intentionally one-helper-per-XML-section so that a future
expansion (parsing angles, dihedrals, nonbonded forces, etc.) can add helpers
alongside the existing ones without disturbing the orchestrator.

## End-to-end workflow context

PDB processing is the second half of the `.off` → OpenMM workflow:

1. `off.openmm.gen_xml(...)` — convert parsed `.off` force-field data into an OpenMM `<ForceField>` XML.
2. `off.openmm.preprocess_pdb(...)` — rename atoms, fix ordering, and add `CONECT` records in the simulation PDB so it matches that XML.

After both steps, the XML and the processed PDB can be loaded directly into an
OpenMM simulation that uses the custom AFM force field.
