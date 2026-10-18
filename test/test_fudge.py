"""Tests for the ``[FUD]`` (1-4 scaling) export: ``off.fudge`` -> ``openmm.gen_xml``.

A pycryoff JSON molecule may carry ``"fudge": {"qq", "vdw", "pairs"}``: the 1-4 scaling its
bonded terms were fitted with. The XML can state one scheme for every 1-4 pair of the system:
``coulomb14scale`` for the charges and the ``afmFudgeQQ`` / ``afmFudgeVDW`` global parameters on
every custom force, which ``afm_openmm.prepare_afm_system`` reads. Anything the XML cannot state
must be refused, never exported as a different model.

Fixture: ``test/sample_json/butanol_intra.json`` (molecule ``UNK``, 1-butanol, ``[Exc]`` = 1-2 +
1-3, i.e. bondCutoff 2), with a card written into a temporary copy.
"""
import copy
import json
import pathlib
import sys
import xml.etree.ElementTree as ET

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
import afmtogmx as afm
from afmtogmx.core import xml_generation as xg

BUTANOL = pathlib.Path(__file__).parent / 'sample_json' / 'butanol_intra.json'


def _one_four(doc, mol):
    off = afm.ReadOFF.from_json(str(BUTANOL))
    return sorted(sorted(p) for p in xg.one_four_pairs(off.bonded, mol))


def _xml(tmp_path, cards, incl_mol=('UNK',)):
    """gen_xml of the butanol fixture with ``cards`` = ``{mol: fudge}`` put into the document."""
    doc = json.loads(BUTANOL.read_text())
    for m in doc['molecules']:
        if m['name'] in cards:
            m['fudge'] = cards[m['name']]
    # The fixture has no per-atom charges; distinct non-zero ones make the charge-product check bite.
    types = sorted({a['cou_type'] for m in doc['molecules'] for a in m['atoms'] if not a.get('pseudo')})
    doc['charges'] = dict(doc.get('charges') or {},
                          per_atom={c: 0.05 * (k + 1) * (-1) ** k for k, c in enumerate(types)})
    src = tmp_path / 'ff.json'
    src.write_text(json.dumps(doc))
    off = afm.ReadOFF.from_json(str(src))
    out = tmp_path / 'ff.xml'
    off.openmm.gen_xml(output=str(out), incl_mol=list(incl_mol))
    return ET.parse(out).getroot(), off


def test_one_four_pairs_of_butanol():
    """1-butanol (C4 O H10): 9 + 9 + 9 + 3 pairs about C-C, C-C, C-C, C-O; none excluded."""
    off = afm.ReadOFF.from_json(str(BUTANOL))
    pairs = xg.one_four_pairs(off.bonded, 'UNK')
    excluded = {frozenset(p) for p in off.bonded['UNK']['EXC']}
    assert pairs and not pairs & excluded
    assert len(pairs) == 30


def test_no_card_writes_the_plain_xml(tmp_path):
    root, _ = _xml(tmp_path, {})
    assert root.find('NonbondedForce').get('coulomb14scale') == '1.0'
    assert not root.findall('.//GlobalParameter')


def test_card_writes_coulomb14scale_and_the_global_parameters(tmp_path):
    card = {'qq': 0.8333, 'vdw': 0.5, 'pairs': _one_four(None, 'UNK')}
    root, off = _xml(tmp_path, {'UNK': card})
    assert off.fudge['UNK']['qq'] == 0.8333
    assert float(root.find('NonbondedForce').get('coulomb14scale')) == 0.8333
    custom = root.findall('CustomNonbondedForce')
    assert custom
    for f in custom:
        g = {p.get('name'): float(p.get('defaultValue')) for p in f.findall('GlobalParameter')}
        assert g == {xg.FUDGE_QQ_PARAM: 0.8333, xg.FUDGE_VDW_PARAM: 0.5}
        assert f.get('bondCutoff') == '2'


def test_card_with_unit_factors_is_a_no_op(tmp_path):
    root, _ = _xml(tmp_path, {'UNK': {'qq': 1.0, 'vdw': 1.0, 'pairs': _one_four(None, 'UNK')}})
    assert root.find('NonbondedForce').get('coulomb14scale') == '1.0'
    assert not root.findall('.//GlobalParameter')


def test_card_listing_a_subset_of_the_1_4_pairs_is_refused(tmp_path):
    card = {'qq': 0.8333, 'vdw': 1.0, 'pairs': _one_four(None, 'UNK')[1:]}
    with pytest.raises(ValueError, match='can only scale every 1-4 pair'):
        _xml(tmp_path, {'UNK': card})


def test_molecule_without_a_card_beside_one_with_a_card_is_refused():
    off = afm.ReadOFF.from_json(str(BUTANOL))
    bonded = copy.deepcopy(off.bonded)
    bonded['B'] = copy.deepcopy(bonded['UNK'])
    fudge = {'UNK': {'qq': 0.8333, 'vdw': 1.0, 'pairs': _one_four(None, 'UNK')}}
    with pytest.raises(ValueError, match='no \\[FUD\\] card'):
        xg.resolve_fudge(bonded, ['UNK', 'B'], fudge, bond_cutoff=2)


def test_molecules_with_different_factors_are_refused():
    off = afm.ReadOFF.from_json(str(BUTANOL))
    bonded = copy.deepcopy(off.bonded)
    bonded['B'] = copy.deepcopy(bonded['UNK'])
    pairs = _one_four(None, 'UNK')
    fudge = {'UNK': {'qq': 0.8333, 'vdw': 1.0, 'pairs': pairs},
             'B': {'qq': 0.5, 'vdw': 1.0, 'pairs': pairs}}
    with pytest.raises(ValueError, match='differ between molecules'):
        xg.resolve_fudge(bonded, ['UNK', 'B'], fudge, bond_cutoff=2)


def test_exclusions_reaching_1_4_are_refused():
    off = afm.ReadOFF.from_json(str(BUTANOL))
    fudge = {'UNK': {'qq': 0.8333, 'vdw': 1.0, 'pairs': _one_four(None, 'UNK')}}
    with pytest.raises(ValueError, match='bondCutoff 3'):
        xg.resolve_fudge(off.bonded, ['UNK'], fudge, bond_cutoff=3)


def test_polarizable_force_field_with_a_card_is_refused():
    off = afm.ReadOFF.from_json(str(BUTANOL))
    fudge = {'UNK': {'qq': 0.8333, 'vdw': 1.0, 'pairs': _one_four(None, 'UNK')}}
    with pytest.raises(NotImplementedError):
        xg.resolve_fudge(off.bonded, ['UNK'], fudge, bond_cutoff=2, polarization={'alphas': {}})


def test_gromacs_writers_refuse_a_card(tmp_path):
    _, off = _xml(tmp_path, {'UNK': {'qq': 0.8333, 'vdw': 1.0, 'pairs': _one_four(None, 'UNK')}})
    with pytest.raises(NotImplementedError, match='GROMACS'):
        off.gmx.gen_bonded_topology(write_to=str(tmp_path / 'b.top'))
    with pytest.raises(NotImplementedError, match='GROMACS'):
        off.gmx.gen_nonbonded_topology(write_to=str(tmp_path / 'nb.top'))


def test_create_system_scales_the_1_4_charge_products(tmp_path):
    """OpenMM's own reading of the XML: every 1-4 exception carries s_qq q_i q_j."""
    openmm = pytest.importorskip('openmm')
    from openmm import app, unit
    pairs = _one_four(None, 'UNK')
    root, off = _xml(tmp_path, {'UNK': {'qq': 0.8333, 'vdw': 1.0, 'pairs': pairs}})
    ff = app.ForceField(str(tmp_path / 'ff.xml'))
    top = app.Topology()
    chain = top.addChain()
    res = top.addResidue('UNK', chain)
    template = ff._templates['UNK']
    atoms = [top.addAtom(a.name, a.element, res) for a in template.atoms]
    for i, j in template.bonds:
        top.addBond(atoms[i], atoms[j])
    system = ff.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    custom = [f for f in system.getForces() if isinstance(f, openmm.CustomNonbondedForce)]
    excluded = {frozenset(custom[0].getExclusionParticles(k)) for k in range(custom[0].getNumExclusions())}
    scaled = 0
    assert any(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) != 0
               for i in range(nb.getNumParticles()))
    for k in range(nb.getNumExceptions()):
        i, j, qq, _, _ = nb.getExceptionParameters(k)
        if frozenset((i, j)) in excluded:
            continue
        qi = nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
        qj = nb.getParticleParameters(j)[0].value_in_unit(unit.elementary_charge)
        assert qq.value_in_unit(unit.elementary_charge ** 2) == pytest.approx(0.8333 * qi * qj, abs=1e-12)
        scaled += 1
    assert scaled == len(pairs)
    for f in custom:
        g = {f.getGlobalParameterName(k): f.getGlobalParameterDefaultValue(k) for k in range(f.getNumGlobalParameters())}
        assert g == {xg.FUDGE_QQ_PARAM: 0.8333, xg.FUDGE_VDW_PARAM: 1.0}
