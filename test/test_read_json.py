"""Tests for the pycryoff JSON reader (``ReadOFF.from_json``, ``core/read_json.py``).

The contract is that a force field loaded from JSON is *indistinguishable* from one parsed out of
a ``.off``: same attribute names, same nesting, same units, so every backend works on it
unchanged. So most of this file checks structural equivalence with the ``.off`` path, and the
rest checks the things JSON carries that a ``.off`` cannot.

Fixtures live in ``test/sample_json/`` and ship with the package:

* ``butanol_intra.json`` — a real pycryoff atomic-force fit of a 1-butanol QM/MM cluster
  (8 bond types, 13 angle types, 4 NCO dihedrals, 77 COU + EXP + SRD nonbonded pairs, two
  molecules). Bonded parameters here are *fitted numbers*, which is the case the ``.off`` path
  gets from its ``Intra-Potential`` section.
* ``extensions_pol.json`` — a synthetic deck exercising what CRYOFF has no syntax for: explicit
  polarizability, a combination rule, an out-of-plane virtual site and declared per-atom charges.
"""
import json
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
import afmtogmx as afm
from afmtogmx.core import functions, read_json

TEST_DIR = pathlib.Path(__file__).parent
SAMPLES = TEST_DIR / 'sample_json'
BUTANOL = SAMPLES / 'butanol_intra.json'
POL = SAMPLES / 'extensions_pol.json'


@pytest.fixture(scope='module')
def butanol():
    return afm.ReadOFF.from_json(str(BUTANOL))


@pytest.fixture(scope='module')
def pol():
    with pytest.warns(UserWarning, match='explicit mutual polarization'):
        return afm.ReadOFF.from_json(str(POL))


# ---------------------------------------------------------------------------------------
# The object is shaped exactly like one built from a .off
# ---------------------------------------------------------------------------------------

def test_bonded_has_the_off_parsers_skeleton(butanol):
    """Every key ``gen_empty_bonded()`` defines must be present, so no consumer KeyErrors."""
    skeleton = functions.gen_empty_bonded()
    for mol, entry in butanol.bonded.items():
        assert set(entry) == set(skeleton), mol
        for card in ('BON', 'ANG', 'BD3', 'DIH', 'CDI'):
            assert set(entry[card]) == set(skeleton[card]), (mol, card)
        assert set(entry['ATO']) == {'All', 'Virtual'}


def test_bonded_is_param_major(butanol):
    """The parameter tuple is the key and the value is every atom group sharing it.

    This inversion is the ``.off`` parser's convention and the GROMACS writers depend on it —
    one ``#define`` per distinct parameter set, then a table of the groups that use it.
    """
    bon = butanol.bonded['UNK']['BON']['HAR']
    assert len(bon) == 8
    for params, groups in bon.items():
        assert isinstance(params, tuple) and len(params) == 2
        assert all(isinstance(v, float) for v in params)
        assert all(len(g) == 2 and all(isinstance(i, int) for i in g) for g in groups)

    # Four dihedral types, one of which covers three methyl-hydrogen groups.
    nco = butanol.bonded['UNK']['DIH']['NCO']
    assert len(nco) == 4
    assert sorted(len(g) for g in nco.values()) == [1, 1, 1, 3]


def test_bond_parameters_are_the_fitted_values(butanol):
    """Native units, straight through: ``re`` in Angstrom, ``kb`` in kcal/mol/A^2."""
    bon = butanol.bonded['UNK']['BON']['HAR']
    by_re = sorted(bon)
    assert by_re[0][0] == pytest.approx(0.9554, abs=1e-3)     # O-H, the shortest bond
    assert by_re[-1][0] == pytest.approx(1.5356, abs=1e-3)    # C-C, the longest
    assert all(200.0 < params[1] < 1500.0 for params in bon)  # plausible stretch constants


def test_angles_and_dihedral_phases_are_in_degrees(butanol):
    """The producer holds radians internally and the schema declares degrees.

    A consistent radians-vs-degrees error round-trips happily and only an absolute check catches
    it — and here it would silently scale every angle force constant's equilibrium by 57.
    """
    assert json.loads(BUTANOL.read_text())['units']['angle'] == 'degree'
    for theta0, _k in butanol.bonded['UNK']['ANG']['HAR']:
        assert 95.0 < theta0 < 125.0, theta0
    for _v, m, delta in butanol.bonded['UNK']['DIH']['NCO']:
        assert m == pytest.approx(3.0)
        assert delta == pytest.approx(0.0)


def test_nonbonded_pairs_are_sorted_and_cou_is_folded(butanol):
    """``.off`` behaviour: the pair key is sorted, and any ``COU*`` label collapses to ``COU``."""
    assert len(butanol.nonbonded) == 77
    for pair, forms in butanol.nonbonded.items():
        assert pair == tuple(sorted(pair))
        assert all(isinstance(sets, list) and sets for sets in forms.values())
    assert butanol.nonbonded[('O0', 'O0')]['COU'] == [[pytest.approx(0.55328307)]]
    assert butanol.nonbonded[('O0', 'O0')]['EXP'][0] == pytest.approx([1075539.9, 4.827992])


def test_exclusions_are_integer_pairs(butanol):
    exc = butanol.bonded['UNK']['EXC']
    assert len(exc) == 39
    assert all(len(p) == 2 and all(isinstance(i, int) for i in p) for p in exc)


def test_backends_and_residues_are_built(butanol):
    """``_finalize`` is shared with the ``.off`` constructor; both paths must end up equipped."""
    assert butanol.gmx is not None and butanol.openmm is not None
    assert set(butanol.residues) == {'Definitions', 'Residues'}
    assert butanol.residues['Definitions']['UNK']['All'][0] == 'O0'


def test_netf_and_torq_are_excluded_everywhere(butanol):
    """The pseudo-atoms are fitting targets, not particles.

    Every filter in this package compares against upper-case ``'NETF'``, which works on a ``.off``
    only because CRYOFF's echo upper-cases them. A deck may write ``NetF``/``Torq``, and the JSON
    preserves that spelling — so the reader normalizes. Without it, two non-atoms reach the charge
    table and the residue definitions.
    """
    assert 'NetF' not in butanol.charges['UNK'] and 'NETF' not in butanol.charges['UNK']
    assert len(butanol.charges['UNK']) == 7
    assert butanol.residues['Definitions']['UNK']['All'] == [
        'O0', 'C1', 'C2', 'C2', 'C3', 'H0', 'H1', 'H1', 'H2', 'H2', 'H2', 'H2', 'H2', 'H2', 'H2']
    assert butanol.residues['Residues']['UNK']['All'] == [list(range(1, 16))]

    types = {t for pair in butanol.nonbonded for t in pair}
    assert not types & {'NETF', 'TORQ', 'NetF', 'Torq'}


# ---------------------------------------------------------------------------------------
# What JSON carries that a .off cannot
# ---------------------------------------------------------------------------------------

def test_declared_charges_arrive_populated(pol):
    """The defect this whole path exists to fix.

    A ``.off`` holds charge *products* only, so ``ReadOFF`` sets every charge to 0.0 and waits for
    a separate ``.charges`` file. Skip that step and the topology is silently neutral everywhere.
    """
    assert pol.charges['MOH'] == {
        'O0': pytest.approx(-0.62089), 'C1': pytest.approx(0.11104),
        'H0': pytest.approx(0.41169), 'H1': pytest.approx(0.04885),
        'LP': pytest.approx(0.0)}


def test_missing_per_atom_charges_warn_rather_than_defaulting_silently():
    with pytest.warns(UserWarning, match='no per-atom charges'):
        off = afm.ReadOFF.from_json(str(BUTANOL))
    assert set(off.charges['UNK'].values()) == {0.0}


def test_out_of_plane_virtual_site_survives(pol):
    """A weighted average is confined to its parents' plane; a lone pair is not."""
    virtual = pol.bonded['MOH']['ATO']['Virtual']
    assert list(virtual) == [(6, 'LP', 'LP')]
    assert virtual[(6, 'LP', 'LP')] == {'kind': 'oop',
                                        'rule': [(0.30, 1), (0.30, 2), (0.55, 3)]}


def test_polarization_is_surfaced_and_warned_about(pol):
    """A model fitted with induction is not the same model without it, so it cannot be dropped."""
    assert pol.polarization['mode'] == 'split'
    assert pol.polarization['poltype'] == 'direct'
    assert pol.polarization['thole'] == pytest.approx(0.42)
    assert pol.polarization['alphas']['O0'] == pytest.approx(0.8552)


def test_combination_rules_are_surfaced(pol):
    assert [(c['target'], c['param'], c['rule']) for c in pol.combinations] == \
           [('EXP', 'alpha', 'harmonic')]


def test_provenance_and_fit_are_available(butanol):
    assert butanol.provenance['ff_file'].endswith('intra.ff')
    assert butanol.fit['kind'] == 'cluster_intra'
    assert butanol.fit['rank'] == 46


# ---------------------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------------------

def test_unknown_schema_is_refused(tmp_path):
    """Best-effort parsing of an unknown schema yields a silently wrong force field."""
    doc = json.loads(POL.read_text())
    doc['schema'] = 'pycryoff-ff/2'
    path = tmp_path / 'future.json'
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match='unsupported force-field schema'):
        afm.ReadOFF.from_json(str(path))


def test_wrong_units_are_refused(tmp_path):
    """kJ/nm instead of kcal/A is wrong by ~4184 with no other symptom."""
    doc = json.loads(POL.read_text())
    doc['units']['energy'] = 'kJ/mol'
    path = tmp_path / 'sikilo.json'
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match='native units'):
        afm.ReadOFF.from_json(str(path))


# ---------------------------------------------------------------------------------------
# The whole point: the backends work untouched
# ---------------------------------------------------------------------------------------

def test_openmm_xml_generation_works(butanol, tmp_path):
    out = tmp_path / 'forcefield.xml'
    butanol.openmm.gen_xml(output=str(out))
    text = out.read_text()
    assert '<ForceField>' in text and '<AtomTypes>' in text
    assert 'UNK_O0' in text and 'CustomBondForce' in text


def test_text_report_works(butanol, tmp_path):
    out = tmp_path / 'report.txt'
    butanol.write_report(str(out))
    text = out.read_text()
    assert 'bonded parameters of UNK' in text
    assert 'U=kb/2*(r-re)^2' in text


def test_off_parsing_is_unaffected():
    """``__init__`` was refactored to share ``_finalize``; the ``.off`` path must be unchanged."""
    off = afm.ReadOFF(str(TEST_DIR / 'sample_off_files' / 'water_intra.off'))
    assert off.charges and off.residues and off.gmx and off.openmm
    assert set(off.bonded['H20QM']) == set(functions.gen_empty_bonded())


# --------------------------------------------------------------------------
# Split vdW / Coulomb atom typing
# --------------------------------------------------------------------------
# A CRYOFF atom line is `<label> <VDWtype> [COUtype]`. When the third column is given, the
# nonbonded cards stay keyed on the vdW type while every OpenMM <Type> comes from the Coulomb
# type -- one repulsion type can serve several charge types. Matching the two namespaces by
# name then silently finds nothing, and because the Discrete2D tables are zero-filled the XML
# comes out with A = 0 for those pairs rather than raising.

_SPLIT_NONBONDED = {('C2', 'C2'): {'EXP': [[30937.8, 3.471]]},
                    ('C2', 'H1'): {'EXP': [[5701.5, 3.394]]},
                    ('O0', 'O0'): {'EXP': [[33520.4, 3.521]]}}
#: five ring carbons sharing one repulsion type, split ortho/meta/para for charges
_SPLIT_BONDED = {'UNK': {'ATO': {'All': {1: ('O0', 'O0'), 2: ('C2', 'C2o'), 3: ('C2', 'C2m'),
                                         4: ('C2', 'C2p'), 5: ('H1', 'H1o'), 6: ('H1', 'H1m')},
                                 'Virtual': {}}}}


def _split_atom_types():
    from afmtogmx.core import xml_generation
    return xml_generation.collect_atom_types(_SPLIT_BONDED, ['UNK'])


def test_atom_types_come_from_the_coulomb_column():
    """OpenMM types must follow the Coulomb type -- that is what carries a per-atom charge."""
    assert [q for _, _, q in _split_atom_types()] == [
        'UNK_O0', 'UNK_C2o', 'UNK_C2m', 'UNK_C2p', 'UNK_H1o', 'UNK_H1m']


def test_nonbonded_expands_across_types_sharing_a_vdw_type():
    """One vdW pair must reach every Coulomb-type pair that maps onto it, same parameters."""
    from afmtogmx.core import xml_generation
    at = _split_atom_types()
    exp, _, _, _ = xml_generation.collect_nonbonded(_SPLIT_NONBONDED, at, _SPLIT_BONDED)

    carbons, hydrogens = ('UNK_C2o', 'UNK_C2m', 'UNK_C2p'), ('UNK_H1o', 'UNK_H1m')
    got = {(q1, q2) for q1, q2, _, _ in exp}
    assert {(a, b) for a in carbons for b in carbons} <= got, 'C2 x C2 not fully expanded'
    assert {(a, b) for a in carbons for b in hydrogens} <= got, 'C2 x H1 not fully expanded'
    assert ('UNK_O0', 'UNK_O0') in got

    shared = {(A, alpha) for q1, q2, A, alpha in exp if q1 in carbons and q2 in carbons}
    assert shared == {(30937.8, 3.471)}, 'expanded pairs must share the fitted parameters'


def test_charges_follow_the_coulomb_column():
    """A charge belongs to the Coulomb type, including when that name is not a vdW type.

    `build_type_to_charge` used to join through the vdW column, which is a no-op on any deck
    whose two type columns are equal and therefore went unnoticed. On a split-typed deck it
    returned 0.0 for every Coulomb type that is not also the name of a vdW type -- here that is
    C2o/C2m/C2p/H1o/H1m, i.e. everything but O0 -- and the emitted XML was complete, well-formed
    and had no electrostatics on most of the molecule. There is no caller-side workaround: three
    charges cannot be stored under the one vdW key `C2`.
    """
    from afmtogmx.core import xml_generation
    at = _split_atom_types()
    charges = {'UNK': {'O0': -0.61, 'C2o': -0.34, 'C2m': -0.06, 'C2p': -0.22,
                       'H1o': 0.16, 'H1m': 0.13}}
    got = xml_generation.build_type_to_charge(_SPLIT_BONDED, charges, at)
    assert got == {'UNK_O0': -0.61, 'UNK_C2o': -0.34, 'UNK_C2m': -0.06,
                   'UNK_C2p': -0.22, 'UNK_H1o': 0.16, 'UNK_H1m': 0.13}
    assert 0.0 not in got.values(), 'a Coulomb type silently lost its charge'


def test_unmatched_vdw_names_warn_instead_of_silently_zeroing():
    """Without `bonded` the pairs are dropped; that must not be silent."""
    from afmtogmx.core import xml_generation
    at = _split_atom_types()
    with pytest.warns(UserWarning, match='match no atom type'):
        exp, _, _, _ = xml_generation.collect_nonbonded(_SPLIT_NONBONDED, at)
    assert [q for q in exp if q[0].startswith('UNK_C2')] == [], 'expected the dropped-pair case'


def test_cou_only_pairs_do_not_warn():
    """Butanol's MM shell interacts by charge alone; an unmatched name there is normal."""
    from afmtogmx.core import xml_generation
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter('error')                       # any warning fails the test
        xml_generation.collect_nonbonded(
            {('O0', 'O01MM'): {'COU': [[-0.0166]]}}, _split_atom_types(), _SPLIT_BONDED)


def test_unsplit_force_fields_are_unaffected():
    """The common case -- no Coulomb column -- must produce exactly what it did before."""
    from afmtogmx.core import xml_generation
    bonded = {'UNK': {'ATO': {'All': {1: ('O0', 'O0'), 2: ('C2', 'C2'), 3: ('H1', 'H1')},
                              'Virtual': {}}}}
    at = xml_generation.collect_atom_types(bonded, ['UNK'])
    without = xml_generation.collect_nonbonded(_SPLIT_NONBONDED, at)[0]
    with_ = xml_generation.collect_nonbonded(_SPLIT_NONBONDED, at, bonded)[0]
    assert sorted(without) == sorted(with_)


def test_butanol_xml_repulsion_table_is_not_silently_zero(butanol, tmp_path):
    """End to end on the real fixture: the QM molecule's repulsion must reach the XML table.

    The MM shell's block is legitimately zero -- it interacts by charge alone -- so this
    checks the QM block specifically rather than the table as a whole.
    """
    import xml.etree.ElementTree as ET
    out = tmp_path / 'forcefield.xml'
    butanol.openmm.gen_xml(output=str(out))
    root = ET.parse(out).getroot()

    types = [t.get('name') for t in root.iter('Type')]
    table = next(f for f in root.iter('Function') if f.get('name') == 'aTable')
    n = int(table.get('xsize'))
    values = [float(v) for v in table.text.split()]
    assert len(values) == n * n == len(types) ** 2

    qm = [i for i, t in enumerate(types) if t.startswith('UNK_')]
    assert qm, 'no QM types in the XML'
    for i in qm:
        for j in qm:
            assert values[i * n + j] != 0.0, f'zero repulsion for {types[i]} x {types[j]}'


# --------------------------------------------------------------------------
# Virtual sites survive the JSON -> OpenMM XML path
# --------------------------------------------------------------------------

def test_json_virtual_site_reaches_the_openmm_xml(tmp_path, monkeypatch):
    """Regression: `from_json` stores the vsite rule structured, the .off parser stores it
    as raw tokens, and the XML writer understood only the tokens -- so a JSON-loaded model
    lost every <VirtualSite> line silently and OpenMM built the site as a free particle."""
    from afmtogmx.core.xml_generation import gen_residues
    bonded = {'W': {'ATO': {'All': {1: ('OW', 'OW'), 2: ('HW', 'HW'),
                                    3: ('HW', 'HW'), 4: ('EW', 'EW')},
                            'Virtual': {(4, 'EW', 'EW'): {'kind': 'average',
                                                          'rule': [(0.6, 1), (0.2, 2), (0.2, 3)]}}},
                    'BON': {}}}
    xml = gen_residues(bonded, ['W'], {})
    assert '<VirtualSite type="average3"' in xml
    assert 'weight1="0.6"' in xml and 'weight2="0.2"' in xml and 'weight3="0.2"' in xml


def test_off_token_virtual_site_still_works():
    """The .off text parser's raw-token form must keep working unchanged."""
    from afmtogmx.core.xml_generation import gen_residues
    bonded = {'W': {'ATO': {'All': {1: ('OW', 'OW'), 2: ('HW', 'HW'),
                                    3: ('HW', 'HW'), 4: ('EW', 'EW')},
                            'Virtual': {(4, 'EW', 'EW'):
                                        ('3:', '0.6', '1', '+', '0.2', '2', '+', '0.2', '3')}},
                    'BON': {}}}
    assert '<VirtualSite type="average3"' in gen_residues(bonded, ['W'], {})


def test_unwritable_virtual_site_raises_instead_of_vanishing():
    """An out-of-plane site is not a linear average; writing average3 would be wrong
    geometry that looks well-formed. Refuse loudly rather than drop the line."""
    from afmtogmx.core.xml_generation import gen_residues
    bonded = {'W': {'ATO': {'All': {1: ('O0', 'O0'), 2: ('C1', 'C1'),
                                    3: ('H0', 'H0'), 6: ('LP', 'LP')},
                            'Virtual': {(6, 'LP', 'LP'): {'kind': 'oop',
                                                          'rule': [(0.3, 1), (0.3, 2), (0.55, 3)]}}},
                    'BON': {}}}
    with pytest.raises(ValueError, match='virtual site'):
        gen_residues(bonded, ['W'], {})


# --------------------------------------------------------------------------
# Charge penetration (CPN) reaches the OpenMM XML with the right energy
# --------------------------------------------------------------------------

FCOV_COU = 332.0637


def _cpn_reference(params, r_angstrom):
    """pycryoff's own CPN form, kcal/mol, r in Angstrom. Transcribed from
    pycryoff.model.functional_forms and used here as the thing to match."""
    import math
    scale, ziqj, zjqi, bi, bj, v0, v1, v2, v3, w0, w1 = params
    r = r_angstrom
    return -FCOV_COU * scale / r * (
        ziqj * (1 + bj * r / 2) * math.exp(-bj * r)
        + zjqi * (1 + bi * r / 2) * math.exp(-bi * r)
        + (v0 + v1 * r + v2 * r ** 2 + v3 * r ** 3) * math.exp(-bi * r)
        + (w0 + w1 * r) * math.exp(-bj * r))


def _energy_from_custom_nb_xml(section, n_types, t1, t2, r_nm):
    """Build the emitted <CustomNonbondedForce> in OpenMM and evaluate it for one pair.

    The point is to run the *string* the generator wrote, not a Python restatement of it.
    """
    import xml.etree.ElementTree as ET
    import openmm as mm
    from openmm import unit

    el = ET.fromstring(section)
    force = mm.CustomNonbondedForce(el.attrib['energy'])
    for fn in el.findall('Function'):
        values = [float(v) for v in fn.text.split()]
        force.addTabulatedFunction(fn.attrib['name'], mm.Discrete2DFunction(
            int(fn.attrib['xsize']), int(fn.attrib['ysize']), values))
    force.addPerParticleParameter('t')
    force.addParticle([t1])
    force.addParticle([t2])

    system = mm.System()
    system.addParticle(1.0)
    system.addParticle(1.0)
    system.addForce(force)
    context = mm.Context(system, mm.VerletIntegrator(1.0),
                         mm.Platform.getPlatformByName('Reference'))
    context.setPositions([(0, 0, 0), (r_nm, 0, 0)])
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)


def test_cpn_xml_reproduces_the_pycryoff_energy_both_ways_round():
    """The CPN form is asymmetric in i and j, and a CustomNonbondedForce picks the order
    of a pair itself. Filling the tables symmetrically would evaluate half the unlike pairs
    with the two valence widths swapped -- a wrong number that looks perfectly ordinary."""
    from afmtogmx.core import xml_generation

    # C0~N0 from an acetonitrile fit: unequal widths, so V and W are both live and the
    # asymmetry is real.
    params = [1.0, -25.20978571, -25.78640356, 3.650181691, 4.12556543,
              -3667.744861, 1078.114446, 0.0, 0.0, 3695.607315, 746.7201562]
    atom_types = [('UNK', 'C0', 'UNK_C0'), ('UNK', 'N0', 'UNK_N0')]
    section = xml_generation.gen_cpn_force([('UNK_C0', 'UNK_N0', params)], atom_types)

    for r_ang in (2.5, 3.0, 3.5, 4.0, 6.0):
        expected = _cpn_reference(params, r_ang) * 4.184      # kcal/mol -> kJ/mol
        forward = _energy_from_custom_nb_xml(section, 2, 0, 1, r_ang * 0.1)
        reverse = _energy_from_custom_nb_xml(section, 2, 1, 0, r_ang * 0.1)
        assert abs(forward - expected) < 1e-6 * max(1.0, abs(expected)), (r_ang, forward, expected)
        assert abs(reverse - expected) < 1e-6 * max(1.0, abs(expected)), (r_ang, reverse, expected)


def test_cpn_like_pair_matches_too():
    """The like-type case (equal widths, cubic V polynomial, W zero) is the other branch of
    pycryoff's closed form and shares no code with the unequal one."""
    from afmtogmx.core import xml_generation

    params = [1.0, -21.02561168, -21.02561168, 3.650181691, 3.650181691,
              23.23800538, 58.31577246, 58.05359043, 23.5451281, 0.0, 0.0]
    atom_types = [('UNK', 'C0', 'UNK_C0')]
    section = xml_generation.gen_cpn_force([('UNK_C0', 'UNK_C0', params)], atom_types)
    for r_ang in (3.0, 4.0, 5.0):
        expected = _cpn_reference(params, r_ang) * 4.184
        got = _energy_from_custom_nb_xml(section, 1, 0, 0, r_ang * 0.1)
        assert abs(got - expected) < 1e-6 * max(1.0, abs(expected)), (r_ang, got, expected)


def test_cpn_vanishes_at_long_range():
    """CPN is a correction that must die away, or it is silently changing the charges."""
    from afmtogmx.core import xml_generation
    params = [1.0, -25.20978571, -25.78640356, 3.650181691, 4.12556543,
              -3667.744861, 1078.114446, 0.0, 0.0, 3695.607315, 746.7201562]
    section = xml_generation.gen_cpn_force(
        [('UNK_C0', 'UNK_N0', params)],
        [('UNK', 'C0', 'UNK_C0'), ('UNK', 'N0', 'UNK_N0')])
    assert abs(_energy_from_custom_nb_xml(section, 2, 0, 1, 1.2)) < 1e-6


# --------------------------------------------------------------------------
# Explicit polarization ([POL]) reaches the OpenMM XML
# --------------------------------------------------------------------------

_POL_BONDED = {
    'UNK': {'ATO': {'All': {1: ('C0', 'C0'), 2: ('H0', 'H0'), 3: ('H0', 'H0'),
                            4: ('H0', 'H0'), 5: ('C1', 'C1'), 6: ('N0', 'N0'),
                            7: ('NETF', 'NETF'), 8: ('TORQ', 'TORQ')},
                    'Virtual': {}},
            'BON': {'HAR': {(1.47, 700.0): [[1, 5]],
                            (1.09, 700.0): [[1, 2], [1, 3], [1, 4]],
                            (1.16, 700.0): [[5, 6]]}},
            'EXC': [[a, b] for a in range(1, 7) for b in range(a + 1, 7)]},
}
_POL_SPEC = {'mode': 'fixed', 'poltype': 'mutual', 'thole': 0.39, 'eps': 1e-05,
             'intramolecular': 'exclude',
             'alphas': {'C0': 1.476114, 'H0': 0.54884, 'C1': 1.476114, 'N0': 1.187309}}


def _pol_atom_types():
    from afmtogmx.core import xml_generation
    return xml_generation.collect_atom_types(_POL_BONDED, ['UNK'])


def test_polarize_types_are_numeric_because_amoeba_parses_them_with_int():
    """OpenMM reads AMOEBA types with int(atom.attrib['type']) and assigns them with
    int(data.atomType[atom]). A descriptive name anywhere in the file is a ValueError at
    createSystem, so the whole file has to be numbered once this force appears."""
    from afmtogmx.core import xml_generation
    at = _pol_atom_types()
    names = xml_generation.build_type_names(at, numeric=True)
    assert sorted(names.values(), key=int) == ['1', '2', '3', '4']
    for name in names.values():
        int(name)
    # the descriptive name survives as the class, which is what every other section uses
    types_xml = xml_generation.gen_atomtypes(_POL_BONDED, at, names)
    assert '<Type name="1" class="UNK_C0"' in types_xml
    assert '<Atom name="C0" type="1"/>' in xml_generation.gen_residues(
        _POL_BONDED, ['UNK'], {}, names)


def test_multipole_force_carries_charges_alphas_and_polarization_groups():
    from afmtogmx.core import xml_generation
    at = _pol_atom_types()
    names = xml_generation.build_type_names(at, numeric=True)
    charges = {'UNK_C0': -0.4589474789, 'UNK_H0': 0.1729011118,
               'UNK_C1': 0.3709080284, 'UNK_N0': -0.430663885}
    xml = xml_generation.gen_multipole_force(_POL_BONDED, ['UNK'], at, names, charges,
                                             _POL_SPEC)
    assert '<Multipole type="1" kz="0" kx="0" c0="-0.4589474789"' in xml
    # Angstrom^3 -> nm^3
    assert 'polarizability="0.001476114"' in xml
    # C0 is bonded to H0 and C1, so both are in its polarization group; the group then
    # closes over the whole molecule, which is what stops it polarizing itself.
    c0_line = [l for l in xml.splitlines() if l.startswith('<Polarize type="1"')][0]
    assert 'pgrp1="2"' in c0_line and 'pgrp2="3"' in c0_line, c0_line
    n0_line = [l for l in xml.splitlines() if l.startswith('<Polarize type="4"')][0]
    assert 'pgrp1="3"' in n0_line
    # every type needs a Polarize tag: the builder's polarizationParams is a defaultdict(dict)
    # and a missing one surfaces as AttributeError deep inside createSystem
    assert xml.count('<Polarize ') == 4


def test_multipole_force_refuses_a_deck_openmm_cannot_express():
    import pytest
    from afmtogmx.core import xml_generation
    at = _pol_atom_types()
    names = xml_generation.build_type_names(at, numeric=True)
    for key, value in (('poltype', 'direct'), ('eps', 1e-8)):
        spec = dict(_POL_SPEC, **{key: value})
        with pytest.raises(ValueError):
            xml_generation.gen_multipole_force(_POL_BONDED, ['UNK'], at, names, {}, spec)


def test_nonbonded_charges_are_zeroed_when_the_multipole_force_has_them():
    """Both forces carrying the charges would double every Coulomb interaction."""
    from afmtogmx.core import xml_generation
    at = _pol_atom_types()
    charges = {'UNK_C0': -0.46, 'UNK_H0': 0.17, 'UNK_C1': 0.37, 'UNK_N0': -0.43}
    xml = xml_generation.gen_nonbonded_force(at, charges, charges_elsewhere=True)
    assert xml.count('charge="0.0"') == 4
    assert '-0.46' not in xml


# --------------------------------------------------------------------------
# bondCutoff is derived from the deck, not assumed
# --------------------------------------------------------------------------

def test_bond_cutoff_follows_the_exclusion_card():
    from afmtogmx.core import xml_generation
    # all 15 pairs excluded -> acetonitrile's widest separation is the 1-4 H...N, so 3
    assert xml_generation.required_bond_cutoff(_POL_BONDED, ['UNK']) == 3

    twelve = {'UNK': dict(_POL_BONDED['UNK'],
                          EXC=[p for p in _POL_BONDED['UNK']['EXC']
                               if tuple(p) not in ((2, 6), (3, 6), (4, 6))])}
    assert xml_generation.required_bond_cutoff(twelve, ['UNK']) == 2


def test_bond_cutoff_refuses_an_exclusion_set_it_cannot_express():
    """An OpenMM CustomNonbondedForce in a ForceField XML can only exclude by bond count.
    Silently using 2 for a deck that means something else exports a different model."""
    import pytest
    from afmtogmx.core import xml_generation
    odd = {'UNK': dict(_POL_BONDED['UNK'], EXC=[[1, 2], [5, 6]])}
    with pytest.raises(ValueError):
        xml_generation.required_bond_cutoff(odd, ['UNK'])
