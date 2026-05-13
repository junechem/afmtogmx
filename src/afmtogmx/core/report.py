"""Publication-ready text report exporter for `ReadOFF`.

Renders the parsed force-field data on a `ReadOFF` instance as a fixed-width
plain-text document suitable for pasting into a Word document at default
margins with a monospace font (Courier New 10pt). Units are the .off file's
native kcal / Angstrom — no conversion is performed.

The public entry point is `write_report(off, path, ...)`. A convenience
method `ReadOFF.write_report(path, ...)` delegates here.
"""

LINE_WIDTH = 95

# ---------------------------------------------------------------------------
# Bonded-subtype column schemas: (header_label, format_spec)
# ---------------------------------------------------------------------------

BONDED_SCHEMAS = {
    'BON': {
        'HAR': {
            'natoms': 2,
            'formula': 'U=kb/2*(r-re)^2',
            'cols': [('re(A)', '>10.4f'), ('kb(kcal/(mol A^2))', '>22.4f')],
        },
        'QUA': {
            'natoms': 2,
            'formula': 'U=quartic(r)',
            'cols': [('P1', '>14.4f'), ('P2', '>14.4f'),
                     ('P3', '>14.4f'), ('P4', '>14.4f')],
        },
    },
    'ANG': {
        'HAR': {
            'natoms': 3,
            'formula': 'U=k0/2*(theta-theta_e)^2',
            'cols': [('theta_e(deg)', '>14.4f'),
                     ('k0(kcal/(mol rad^2))', '>22.4f')],
        },
        'QUA': {
            'natoms': 3,
            'formula': 'U=quartic(theta)',
            'cols': [('P1', '>14.4f'), ('P2', '>14.4f'),
                     ('P3', '>14.4f'), ('P4', '>14.4f')],
        },
    },
    'BD3': {
        'QBB': {
            'natoms': 3,
            'formula': 'U=QBB out-of-plane',
            'cols': [('P1', '>12.4f'), ('P2', '>12.4f'),
                     ('P3', '>12.4f'), ('P4', '>12.4f'), ('P5', '>12.4f')],
        },
        'MUB': {
            'natoms': 3,
            'formula': 'U=MUB out-of-plane',
            'cols': [('P1', '>12.4f'), ('P2', '>12.4f'),
                     ('P3', '>12.4f'), ('P4', '>12.4f')],
        },
    },
    'DIH': {
        'HAR': {
            'natoms': 4,
            'formula': 'U=k/2*(phi-phi_e)^2',
            'cols': [('P1', '>14.4f'), ('P2', '>14.4f')],
        },
        'NCO': {
            'natoms': 4,
            # Stored tuple is (V_D, m, delta); we render in PDF order: delta, V_D, m.
            'formula': 'U=V_D*(1+cos(m*w-delta))',
            'cols': [('delta(deg)', '>12.4f'),
                     ('V_D(kcal/mol)', '>16.4f'),
                     ('m', '>4.0f')],
            'remap': lambda P: (P[2], P[0], P[1]),
        },
        'COS': {
            'natoms': 4,
            'formula': 'U=V_D*cos(m*w-delta)',
            'cols': [('delta(deg)', '>12.4f'),
                     ('V_D(kcal/mol)', '>16.4f'),
                     ('m', '>4.0f')],
            'remap': lambda P: (P[2], P[0], P[1]),
        },
    },
    # Both 'CDI' and 'CDIH' get the same schema — _parse_bonded actually stores
    # combined dihedrals under 'CDIH', while gen_empty_bonded preallocates 'CDI'.
    'CDI': {
        'CNCO': {
            'natoms': 4,
            'formula': 'U=CNCO combined dihedral',
            'cols': [('P1', '>12.4f'), ('P2', '>12.4f'),
                     ('P3', '>12.4f'), ('P4', '>12.4f')],
        },
        'CCOS': {
            'natoms': 4,
            'formula': 'U=CCOS combined dihedral',
            'cols': [('P1', '>12.4f'), ('P2', '>12.4f'),
                     ('P3', '>12.4f'), ('P4', '>12.4f')],
        },
    },
}
BONDED_SCHEMAS['CDIH'] = BONDED_SCHEMAS['CDI']

BONDED_CATEGORY_LABEL = {
    'BON': 'Bonds',
    'ANG': 'Angles',
    'BD3': 'BD3 out-of-plane',
    'DIH': 'Dihedral',
    'CDI': 'Combined dihedrals',
    'CDIH': 'Combined dihedrals',
}

# Order in which bonded categories appear in the report.
BONDED_CATEGORY_ORDER = ['BON', 'ANG', 'BD3', 'DIH', 'CDI', 'CDIH']

# ---------------------------------------------------------------------------
# Nonbonded interaction layout
# ---------------------------------------------------------------------------
#
# The CRYOFF manual (page 23, Table 6 footnote) is explicit:
#
#   "The CRYOFF symbols are the minimal symbols to be used in the .ff file to
#    denote their functional forms; Addition letters could follow the minimal
#    symbol for identification such as EXPinter."
#
# So a stored key like 'EXPINTRA' / 'EXPINTER' / 'EXPW' / 'STRC' / 'STRCINTER'
# is the user-chosen full identifier of a fit, and its first 3 characters
# (EXP / STR / ...) name the functional form. The report groups every stored
# key under its 3-char-prefix section.
#
# Two label conventions are supported:
#   notation='standard' (default): publication-style names from the manual —
#                                  A, alpha, Cn, R0 — consistent within a row.
#   notation='PN':                 generic P1, P2, P3, ... matching the
#                                  manual's positional parameter order.
#
# For POW and SRD the exponent stored by CRYOFF is signed (e.g. -6 for
# attractive dispersion). The 'standard' notation reports the magnitude in
# the `n` column and folds that magnitude into the coefficient label as
# `C{n}` (e.g. C6), matching the convention used in Wang-group publications.
# The coefficient itself keeps its sign, so the formula `U = Cn / r^n` with
# positive n recovers the correct attractive/repulsive sign.

# Section render order: Coulombics first, then short-range repulsion/dispersion
# forms following the manual's Table 6 layout.
NONBONDED_ORDER = ['COU', 'THC',
                   'GLJ', 'BUC', 'DBU', 'STR', 'EXP', 'GEX',
                   'POW', 'CSP', 'GDP', 'PEX', 'DPO', 'SRD', 'TTP']


def _prefix(itype):
    """Return the 3-char minimal symbol for a stored interaction-type key.

    Per the CRYOFF manual (Table 6 footnote, page 23), the first 3 characters
    of an interaction-type name determine the functional form; any trailing
    characters are a user-chosen identifier (e.g. EXPINTRA, EXPINTER, EXPW,
    STRC, STRCINTER). All such suffixed names share the prefix's schema.
    """
    return itype[:3] if isinstance(itype, str) else itype


def _nonbonded_layout(prefix, rows, notation):
    """Return ``(formula, [(header, fmt), ...], row_transform)`` for a section.

    The row_transform takes a stored param list and returns the values to feed
    into the column format specs. For most prefixes it is the identity. For
    POW/SRD/STR under standard notation it folds the signed exponent in the
    stored data into the column's magnitude convention.
    """
    if notation == 'PN':
        return _pn_layout(prefix, rows)

    if prefix == 'COU':
        return ('COU = q1*q2 / (4*pi*eps0*r)',
                [('q1*q2(e^2)', '>14.5f')],
                lambda r: r)
    if prefix == 'THC':
        return ('THC = (q1*q2 / (4*pi*eps0*r)) * Thole_damping(alpha)',
                [('q1*q2(e^2)', '>14.5f'), ('alpha', '>12.4f')],
                lambda r: r)
    if prefix == 'EXP':
        return ('EXP = A * exp(-alpha * r)',
                [('A(kcal/mol)', '>16.3f'), ('alpha(1/A)', '>14.3f')],
                lambda r: r)
    if prefix == 'BUC':
        return ('BUC = A * exp(-alpha * r) + C6 / r^6',
                [('A(kcal/mol)', '>16.3f'),
                 ('C6(kcal mol^-1 A^6)', '>22.3f'),
                 ('alpha(1/A)', '>14.3f')],
                lambda r: r)
    if prefix == 'GLJ':
        return ('GLJ = C(n1)/r^n1 + C(n2)/r^n2',
                [('C(n1)(kcal/mol)', '>16.4f'), ('C(n2)(kcal/mol)', '>16.4f'),
                 ('n1', '>5.2f'), ('n2', '>5.2f')],
                lambda r: r)
    if prefix == 'POW':
        return ('POW = Cn / r^n',
                [('Cn(kcal/mol)', '>16.3f'), ('n', '>4.0f')],
                lambda r: (r[0], abs(r[1])))
    if prefix == 'SRD':
        return ('SRD = Cn / (r^n + R0^n)',
                [('Cn(kcal/mol)', '>16.3f'), ('n', '>4.0f'), ('R0(A)', '>10.3f')],
                lambda r: (r[0], abs(r[1]), r[2]))
    if prefix == 'STR':
        formula = 'STR = A * [(1/r^n - 1/R0^n) + n*(r-R0)/R0^(n+1)]  for r<=R0  (shift-truncated)'
        return (formula,
                [('A(kcal/mol)', '>16.3f'), ('n', '>4.0f'), ('R0(A)', '>10.3f')],
                lambda r: (r[0], abs(r[1]), r[2]))
    if prefix == 'PEX':
        return ('PEX = A * r^n * exp(-alpha * r)',
                [('A(kcal/mol)', '>16.3f'), ('n', '>5.2f'),
                 ('alpha(1/A)', '>14.3f')],
                lambda r: (r[0], abs(r[1]), r[2]) if len(r) >= 3 else r)
    if prefix == 'GDP':
        return ('GDP = A * r^n * exp(-beta * r^2)',
                [('A(kcal/mol)', '>16.3f'), ('n', '>5.2f'),
                 ('beta(1/A^2)', '>14.4f')],
                lambda r: (r[0], abs(r[1]), r[2]) if len(r) >= 3 else r)
    if prefix == 'DPO':
        return ('DPO = A * r^n / (1 + exp(-beta * (r - r0)))',
                [('A(kcal/mol)', '>16.3f'), ('n', '>5.2f'),
                 ('beta(1/A)', '>10.3f'), ('r0(A)', '>10.3f')],
                lambda r: (r[0], abs(r[1]), r[2], r[3]) if len(r) >= 4 else r)
    if prefix == 'GEX':
        return ('GEX = (1 + a*r + b*r^2)^(-c)',
                [('P1', '>14.4f'), ('a', '>12.4f'),
                 ('b', '>12.4f'), ('c', '>8.3f')],
                lambda r: r)
    if prefix == 'CSP':
        return ('CSP = (A/2) * r^n * [cos(pi*(r-r1)/(r2-r1)) + 1]  for r1<r<r2; 0 otherwise',
                [('A(kcal/mol)', '>16.3f'), ('n', '>5.2f'),
                 ('r1(A)', '>10.3f'), ('r2(A)', '>10.3f')],
                lambda r: (r[0], abs(r[1]), r[2], r[3]) if len(r) >= 4 else r)
    if prefix == 'DBU':
        return ('DBU = A * exp(-alpha * r) + C6 * f(r);  f = 1/(1 + exp(-beta*(r-r0)))',
                [('A(kcal/mol)', '>16.3f'),
                 ('C6(kcal mol^-1 A^6)', '>22.3f'),
                 ('alpha(1/A)', '>10.3f'),
                 ('beta(1/A)', '>10.3f'),
                 ('r0(A)', '>10.3f')],
                lambda r: r)
    if prefix == 'TTP':
        return ('TTP = A * r^n * [1 - exp(-beta r) * sum_{k=0..n} (beta r)^k / k!]',
                [('A(kcal/mol)', '>16.3f'), ('n', '>5.2f'),
                 ('beta(1/A)', '>10.3f')],
                lambda r: (r[0], abs(r[1]), r[2]) if len(r) >= 3 else r)
    # Unknown prefix — fall through to the PN layout.
    return _pn_layout(prefix, rows)


def _pn_layout(prefix, rows):
    """Generic PN layout: P1, P2, P3, ... sized to the longest row."""
    max_len = max((len(r) for r in rows), default=0)
    cols = [(f'P{i+1}', '>14.4f') for i in range(max_len)]
    formula = f'{prefix} = function of (P1, P2, ..., P{max_len})' if max_len else f'{prefix} = (custom)'
    return formula, cols, lambda r: r

ATOM_W = 6  # column width for atom-name cells


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def write_report(off, path, incl_mol=None, name_translation=None,
                 molname_translation=None, notation='standard'):
    """Write a publication-style text report of an :class:`ReadOFF` object.

    Parameters
    ----------
    off : ReadOFF
        Source object holding parsed bonded, nonbonded, and charge data.
    path : str
        Destination file path. The file is overwritten if it already exists.
    incl_mol : list of str, optional
        Molecule names to include. ``None`` or an empty list means all.
        Nonbonded pairs are included under a molecule's section if at least
        one of the two atom names belongs to that molecule, so cross-molecule
        (e.g. solute-solvent) interactions are preserved.
    name_translation : dict, optional
        ``{off_atom_name: display_name}``. Applied to every atom name in
        every table. Atoms not in the dict are rendered unchanged.
    molname_translation : dict, optional
        ``{off_molecule_name: display_name}``. Applied to molecule names in
        banners and section headers. Molecules not in the dict are rendered
        unchanged.
    notation : {'standard', 'PN'}, default 'standard'
        Column label convention for nonbonded sections.
        - 'standard': publication-style names (A, alpha, Cn, R0, ...). For
          POW and SRD, the coefficient column is ``C{|n|}`` (e.g. ``C6``)
          and the exponent column reports ``|n|``; the coefficient keeps its
          fitted sign so ``Cn / r^n`` recovers the correct attractive vs.
          repulsive contribution.
        - 'PN': generic ``P1, P2, P3, ...`` matching the manual's positional
          ordering, with no transformation of stored values.
    """
    if notation not in ('standard', 'PN'):
        raise ValueError(
            f"notation must be 'standard' or 'PN', got {notation!r}"
        )
    if not incl_mol:
        incl_mol = list(off.bonded.keys())
    name_translation = dict(name_translation or {})
    molname_translation = dict(molname_translation or {})

    issues = []
    out = []

    for molname in off.bonded:
        if molname not in incl_mol:
            continue
        out.extend(_render_molecule(
            off, molname, name_translation, molname_translation, issues,
            notation,
        ))
        out.append('')

    if issues:
        msg = (
            "write_report refusing to emit: the following items could not "
            "be resolved against off.bonded/off.nonbonded:\n  - "
            + "\n  - ".join(issues)
        )
        raise ValueError(msg)

    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(out))


# ---------------------------------------------------------------------------
# Per-molecule rendering
# ---------------------------------------------------------------------------


def _render_molecule(off, molname, name_tr, mol_tr, issues, notation):
    display = mol_tr.get(molname, molname)
    lines = []

    # Banner
    lines.append('=' * LINE_WIDTH)
    title = f'AFM Force field parameters ({display})'
    lines.append(title.center(LINE_WIDTH))
    lines.append('=' * LINE_WIDTH)

    # Charges
    lines.extend(_render_charges(off.charges.get(molname, {}), name_tr))

    # Bonded
    lines.append('=' * LINE_WIDTH)
    lines.append(f'bonded parameters of {display}:')
    lines.append('=' * LINE_WIDTH)
    lines.append('Parameters are in the .off file native units (kcal/mol, Angstrom).')
    lines.append('')

    bonded = off.bonded.get(molname, {})
    ato_map = bonded.get('ATO', {}).get('All', {})

    rendered_any_bonded = False
    for category in BONDED_CATEGORY_ORDER:
        if category not in bonded:
            continue
        block = _render_bonded_category(
            bonded[category], category, ato_map, molname, name_tr, issues
        )
        if block:
            rendered_any_bonded = True
            lines.extend(block)
    if not rendered_any_bonded:
        lines.append('(no bonded interactions)')
        lines.append('')

    # Nonbonded
    lines.append('=' * LINE_WIDTH)
    lines.append(f'nonbonded parameters of {display}:')
    lines.append('=' * LINE_WIDTH)

    mol_atom_names = set(off.charges.get(molname, {}).keys())
    relevant_pairs = _select_nonbonded_pairs(off.nonbonded, mol_atom_names)
    lines.extend(_render_nonbonded(relevant_pairs, name_tr, notation))

    lines.append('=' * LINE_WIDTH)
    return lines


# ---------------------------------------------------------------------------
# Charges
# ---------------------------------------------------------------------------


def _render_charges(charges_for_mol, name_tr):
    if not charges_for_mol:
        return []
    header = f'{"Type":<{ATOM_W}}  Partial-charges(e)'
    rule = '-' * len(header)
    rows = [header, rule]
    for atom_name, q in charges_for_mol.items():
        display = name_tr.get(atom_name, atom_name)
        rows.append(f'{display:<{ATOM_W}}  {q:>10.5f}')
    return rows


# ---------------------------------------------------------------------------
# Bonded
# ---------------------------------------------------------------------------


def _render_bonded_category(category_dict, category, ato_map, molname,
                            name_tr, issues):
    """Render one category (BON/ANG/BD3/DIH/CDI/CDIH) across its subtypes."""
    schemas = BONDED_SCHEMAS[category]
    label = BONDED_CATEGORY_LABEL[category]

    subtype_blocks = []
    for subtype, params_dict in category_dict.items():
        if not params_dict:
            continue
        if subtype not in schemas:
            issues.append(
                f"Unknown {category} subtype '{subtype}' for molecule "
                f"'{molname}' (no rendering schema)"
            )
            continue
        block = _render_bonded_subtype(
            label, subtype, schemas[subtype], params_dict,
            ato_map, molname, name_tr, issues,
        )
        subtype_blocks.append(block)
    return [line for block in subtype_blocks for line in block]


def _render_bonded_subtype(label, subtype, schema, params_dict, ato_map,
                           molname, name_tr, issues):
    natoms = schema['natoms']
    formula = schema['formula']
    cols = schema['cols']
    remap = schema.get('remap')

    # Build header line: label (subtype)  formula   col1   col2 ...
    atom_header = ''.join(f'{f"at{i+1}":<{ATOM_W}}' for i in range(natoms))
    col_header = ''.join(f'{name:>{_fmt_width(spec)}}' for name, spec in cols)
    header = f'{label} ({subtype})  {formula}'
    columns_line = f'{atom_header}{col_header}'

    lines = [header, columns_line, '-' * max(len(header), len(columns_line))]

    seen = set()
    for param_tuple, atom_lists in params_dict.items():
        rendered_params = remap(param_tuple) if remap else param_tuple
        for atnums in atom_lists:
            try:
                names = [_atom_name(ato_map, n) for n in atnums]
            except KeyError as e:
                issues.append(
                    f"Unknown atom number {e.args[0]} in {category_of(label)} "
                    f"({subtype}) for molecule '{molname}'"
                )
                continue
            sig = _canonical_atom_signature(names)
            dedupe_key = (sig, param_tuple)
            if dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            display_names = [name_tr.get(n, n) for n in names]
            atom_cells = ''.join(f'{n:<{ATOM_W}}' for n in display_names)
            param_cells = ''.join(
                f'{rendered_params[i]:{spec}}' for i, (_, spec) in enumerate(cols)
            )
            lines.append(f'{atom_cells}{param_cells}')

    lines.append('')
    return lines


def _canonical_atom_signature(names):
    """Return a signature equal for reversible orderings of the same term.

    Bonds (2 atoms): order-independent (frozenset).
    Angles (3 atoms): central atom fixed; ends order-independent.
    Dihedrals (4 atoms): equal to its reverse.
    BD3 (3 atoms): treated like angles.
    Other lengths: tuple as-is.
    """
    n = len(names)
    if n == 2:
        return ('bond', tuple(sorted(names)))
    if n == 3:
        return ('angle', names[1], tuple(sorted((names[0], names[2]))))
    if n == 4:
        fwd = tuple(names)
        rev = fwd[::-1]
        return ('dihedral', min(fwd, rev))
    return tuple(names)


def category_of(label):
    """Inverse of BONDED_CATEGORY_LABEL for error messages."""
    for k, v in BONDED_CATEGORY_LABEL.items():
        if v == label:
            return k
    return label


def _atom_name(ato_map, atnum):
    """Resolve an atom number to its Coulomb-type name via ATO['All']."""
    # ATO['All'] keys may be int or str depending on parser path; try both.
    if atnum in ato_map:
        entry = ato_map[atnum]
    elif str(atnum) in ato_map:
        entry = ato_map[str(atnum)]
    else:
        raise KeyError(atnum)
    return entry[1]  # (vdw_type, coulomb_type) → coulomb_type


def _fmt_width(spec):
    """Extract width digits from a format spec like '>16.3f' → 16."""
    digits = ''
    for ch in spec:
        if ch.isdigit():
            digits += ch
        elif digits:
            break
    return int(digits) if digits else 0


# ---------------------------------------------------------------------------
# Nonbonded
# ---------------------------------------------------------------------------


def _select_nonbonded_pairs(nonbonded, mol_atom_names):
    """Return dict of pairs where at least one atom belongs to the molecule."""
    return {
        pair: ints for pair, ints in nonbonded.items()
        if pair[0] in mol_atom_names or pair[1] in mol_atom_names
    }


def _render_nonbonded(pairs, name_tr, notation):
    if not pairs:
        return ['(no nonbonded interactions)', '']

    # Discover which 3-char prefixes appear across all stored keys. A stored
    # key like 'EXPINTRA' contributes prefix 'EXP'; 'STRC' contributes 'STR'.
    # Preserve first-appearance order for any prefix not in NONBONDED_ORDER.
    present_prefixes = []
    seen = set()
    rows_by_prefix = {}
    for ints in pairs.values():
        for itype, plist in ints.items():
            pref = _prefix(itype)
            if pref not in seen:
                seen.add(pref)
                present_prefixes.append(pref)
            rows_by_prefix.setdefault(pref, []).extend(plist)

    ordered = [p for p in NONBONDED_ORDER if p in seen]
    extras = [p for p in present_prefixes if p not in ordered]
    ordered.extend(extras)

    # Resolve each section's layout once so the energy-expressions block at
    # the top of the nonbonded summary matches the column headers below.
    layouts = {
        pref: _nonbonded_layout(pref, rows_by_prefix.get(pref, []), notation)
        for pref in ordered
    }

    lines = ['Energy expressions:']
    for pref in ordered:
        formula, _, _ = layouts[pref]
        lines.append(f'  {formula}')
    lines.append('')

    for pref in ordered:
        lines.extend(
            _render_nonbonded_section(pref, pairs, name_tr, layouts[pref])
        )
    return lines


def _render_nonbonded_section(prefix, pairs, name_tr, layout):
    """Render every stored-key row whose 3-char prefix matches `prefix`."""
    _, cols, row_transform = layout

    # Formula already shown once in the "Energy expressions:" block above; keep
    # the per-section header short so total line width stays within 95 chars.
    header = f'{prefix} interactions:'
    atom_header = f'{"atom1":<{ATOM_W}}{"atom2":<{ATOM_W}}'
    col_header = ''.join(f'{name:>{_fmt_width(spec)}}' for name, spec in cols)
    columns_line = f'{atom_header}{col_header}'

    block = [header, columns_line, '-' * max(len(header), len(columns_line))]

    any_rows = False
    for (a1, a2), ints in pairs.items():
        for itype, param_list in ints.items():
            if _prefix(itype) != prefix:
                continue
            for params in param_list:
                values = row_transform(params)
                d1 = name_tr.get(a1, a1)
                d2 = name_tr.get(a2, a2)
                atom_cells = f'{d1:<{ATOM_W}}{d2:<{ATOM_W}}'
                param_cells = ''.join(
                    f'{values[i]:{spec}}' if i < len(values)
                    else f'{"":>{_fmt_width(spec)}}'
                    for i, (_, spec) in enumerate(cols)
                )
                block.append(f'{atom_cells}{param_cells}')
                any_rows = True

    if not any_rows:
        block.append('(none for selected molecule)')
    block.append('')
    return block
