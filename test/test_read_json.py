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
