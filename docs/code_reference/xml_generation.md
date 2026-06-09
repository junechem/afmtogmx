# `xml_generation.py` — OpenMM `<ForceField>` XML section builders

Module-level helpers that turn `ReadOFF` data (`bonded`, `nonbonded`, `charges`)
into OpenMM XML sections. Used by `OpenMMBackend.gen_xml`. Atom types are
**qualified** by molecule name (`"UNK_C1"`) so two molecules reusing a raw type
label stay globally unique.

**OpenMM dependency:** tries to import `openmm.app.Element`; if unavailable,
`_OPENMM_AVAILABLE = False` and atom masses fall back to 0.0.

## Internal helpers
- `_atom_map(ato_all)` — normalize `ATO['All']` to `{int_id: (atname, attype)}`.
- `_qualify(mol, attype)` — `"<mol>_<attype>"`.
- `_is_virtual(bonded, mol, attype)` — True if the type is a virtual site in `mol`.
- `_element_symbol(attype)` — element symbol from an atom-type name. Tries the
  first **two** characters as an element symbol first (so `NA`→`Na`, `CL`→`Cl`),
  then the first character (`OW`→`O`); returns the canonical OpenMM symbol on a
  match, else `attype[0]`. Falls back to `attype[0]` when OpenMM is unavailable.
- `_get_mass(element_symbol)` — element mass from OpenMM (0.0 if unavailable).
- `_unique_atom_names(atom_map)` — per-residue unique names `O0, H0, Na0, …`
  (element symbol from `_element_symbol` + counter); skips NETF/TORQ.
- `_matrix(n, val=0.0)` / `_matrix_str(matrix)` — N×N table helpers for the
  Discrete2D lookup tables used by the custom nonbonded forces.

## Data collection (called by the orchestrator)
- `collect_atom_types(bonded, mol_names)` — `[(mol, raw_type, qualified_type), …]`
  in ATO appearance order; skips NETF/TORQ.
- `build_type_to_charge(bonded, charges, atom_types)` — `{qualified_type: charge}`
  by joining ATO atom names to the charges dict.
- `collect_nonbonded(nonbonded, atom_types)` — splits nonbonded entries by type and
  expands raw type names to all matching qualified types. Returns
  `(exp_entries, str_entries, srd_by_power)` in kcal/Å units. POW → SRD with r0=0;
  BUC → one EXP entry + one SRD(power=-6) entry.

## Section builders (each returns a string, or `''` if empty)
- `gen_atomtypes(bonded, atom_types)` — `<AtomTypes>`; element/mass come from
  `_element_symbol` (handles two-letter elements like Na/Cl); virtual sites get
  mass 0.0 and no element.
- `gen_residues(bonded, mol_names, molname_translations)` — `<Residues>` with
  `<Atom>`, `<Bond>` (from `BON` data, deduped), and `<VirtualSite>` entries.
  Residue name comes from `molname_translations` (falls back to molname).
- `_virtual_site_xml(site_name, definition, unique)` — parse a virtual-site tuple
  into an `average2`/`average3` `<VirtualSite>` element (or `None`).
- `gen_nonbonded_force(atom_types, type_to_charge)` — `<NonbondedForce>` carrying
  point charges only (sigma=epsilon=0).
- `gen_bond_force(bonded, mol_names)` — unified `<CustomBondForce>` for HAR
  (k3=k4=0) and QUA bonds, with the quartic energy expression. Converts r0 (×0.1),
  k2 (×418.4), k3 (×4184), k4 (×41840).
- `_bond_pair_qualified(mol, at1_id, at2_id, atom_map)` — `(qual1, qual2)` or
  `(None, None)` if an atom is missing/NETF/TORQ.
- `gen_angle_force(bonded, mol_names)` — `<CustomAngleForce>` for HAR angles
  (theta0→radians, k×4.184).
- `gen_dihedral_force(bonded, mol_names)` — `<PeriodicTorsionForce>` for NCO
  dihedrals (k kcal→kJ; phase already radians; reverse-duplicate suppressed).
- `gen_exp_force(exp_entries, atom_types)` — `<CustomNonbondedForce>` for
  `U=A·exp(-alpha·r)` via Discrete2D tables (A×4.184, alpha×10).
- `gen_srd_force(entries, power, atom_types)` — `<CustomNonbondedForce>` for
  `U=disp/(r^|p|+r0^|p|)` (disp scaled by `4.184·0.1^-power`, r0×0.1).
- `gen_str_force(str_entries, atom_types)` — `<CustomNonbondedForce>` for the
  shifted-truncated power potential; guards against `0^0` for unused pairs.
- `_custom_nb_xml(energy, tables, qualified, n)` — assembles a
  `<CustomNonbondedForce>` block with N×N Discrete2D `<Function>` lookup tables and
  one `<Atom type=… t=index/>` per qualified type.

## Output
- `write_xml(path, sections)` — wraps the section strings in `<ForceField>…</ForceField>`
  and writes to `path`.
