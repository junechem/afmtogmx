"""
Regression tests for afmtogmx.

Each test replicates the corresponding generate_testN.py workflow, writes
outputs to a pytest tmp_path, and compares every output file against the
stored baselines in test/baseline_outputs/testN_*/.

  .xvg files — numeric comparison (numpy.allclose, rtol=1e-6)
  .top files — line-by-line text comparison (trailing whitespace stripped)

A failing test means a code change altered outputs — either a bug was
introduced, or the baselines need to be regenerated intentionally.
"""
import os
import sys
import shutil
import pathlib

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
import afmtogmx as afm

TEST_DIR    = pathlib.Path(__file__).parent
BASELINE    = TEST_DIR / 'baseline_outputs'
SAMPLES     = TEST_DIR / 'sample_off_files'

# Files that are inputs to the workflow, not outputs — skip when comparing
_INPUTS = {'template.top', 'charges.txt'}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def restore_cwd():
    """Restore the working directory after each test."""
    original = os.getcwd()
    yield
    os.chdir(original)


def _setup(tmp_path: pathlib.Path, baseline_name: str, extra_inputs=()) -> pathlib.Path:
    """Copy input files to tmp_path, chdir there, and return the baseline dir."""
    baseline_dir = BASELINE / baseline_name
    shutil.copy(baseline_dir / 'template.top', tmp_path / 'template.top')
    for fname in extra_inputs:
        shutil.copy(baseline_dir / fname, tmp_path / fname)
    os.chdir(tmp_path)
    return baseline_dir


def _compare_xvg(actual: pathlib.Path, baseline: pathlib.Path, rtol=1e-6, atol=1e-10):
    act  = np.loadtxt(actual)
    base = np.loadtxt(baseline)
    assert act.shape == base.shape, (
        f"Shape mismatch in {actual.name}: got {act.shape}, expected {base.shape}"
    )
    if not np.allclose(act, base, rtol=rtol, atol=atol, equal_nan=True):
        diff = np.abs(act - base)
        idx  = np.unravel_index(np.argmax(diff), diff.shape)
        pytest.fail(
            f"Numeric mismatch in {actual.name}\n"
            f"  max abs diff: {diff[idx]:.3e}  at row {idx[0]}, col {idx[1]}\n"
            f"  actual={act[idx]:.6e}  baseline={base[idx]:.6e}"
        )


def _compare_top(actual: pathlib.Path, baseline: pathlib.Path):
    act_lines  = [l.rstrip() for l in actual.read_text().splitlines()]
    base_lines = [l.rstrip() for l in baseline.read_text().splitlines()]
    assert act_lines == base_lines, (
        f"Content mismatch in {actual.name}  "
        f"({len(act_lines)} vs {len(base_lines)} lines)"
    )


def _compare_all(tmp_path: pathlib.Path, baseline_dir: pathlib.Path):
    """Walk baseline_dir; compare every .xvg and output .top against tmp_path."""
    for baseline_file in sorted(baseline_dir.rglob('*')):
        if not baseline_file.is_file():
            continue
        if baseline_file.suffix == '.py':
            continue
        if baseline_file.name in _INPUTS:
            continue

        rel         = baseline_file.relative_to(baseline_dir)
        actual_file = tmp_path / rel

        assert actual_file.exists(), f"Expected output not generated: {rel}"

        if baseline_file.suffix == '.xvg':
            _compare_xvg(actual_file, baseline_file)
        elif baseline_file.suffix == '.top':
            _compare_top(actual_file, baseline_file)


# ---------------------------------------------------------------------------
# Tests 1-8  (test9 is in-progress and excluded)
# ---------------------------------------------------------------------------

def test1_methane_basic(tmp_path):
    """Basic methane workflow — nonbonded tabpot + topology, no special options."""
    baseline = _setup(tmp_path, 'test1_methane_basic')

    off = afm.ReadOFF(off_loc=str(SAMPLES / 'methane_intra.off'))
    nonbonded_tabpot = off.gen_nonbonded_tabpot()
    off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, tabpot_prefix='table')
    off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
    off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')

    _compare_all(tmp_path, baseline)


def test2_ethane_bonded(tmp_path):
    """Ethane workflow — nonbonded + bonded tabpots + topology."""
    baseline = _setup(tmp_path, 'test2_ethane_bonded')

    off = afm.ReadOFF(off_loc=str(SAMPLES / 'ethane_intra.off'))
    nonbonded_tabpot = off.gen_nonbonded_tabpot()
    off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, tabpot_prefix='table')
    bonded_tabpot = off.gen_bonded_tabpot()
    off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot, tabpot_prefix='table')
    off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
    off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')

    _compare_all(tmp_path, baseline)


def test3_water_charges(tmp_path):
    """Water workflow — manual charge assignment, incl_mol filter."""
    baseline = _setup(tmp_path, 'test3_water_charges')

    off = afm.ReadOFF(off_loc=str(SAMPLES / 'water_intra.off'))
    if 'H20QM' in off.charges:
        off.charges['H20QM']['OQM'] = -0.82
        off.charges['H20QM']['HQM'] =  0.41
        off.charges['H20QM']['EQM'] =  0.0

    nonbonded_tabpot = off.gen_nonbonded_tabpot(incl_mol=['H20QM'])
    off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, tabpot_prefix='table')
    bonded_tabpot = off.gen_bonded_tabpot(incl_mol=['H20QM'])
    off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot, tabpot_prefix='table')
    off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top',
                               incl_mol=['H20QM'])
    off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top',
                            incl_mol=['H20QM'], bonded_tabpot=bonded_tabpot)

    _compare_all(tmp_path, baseline)


def test4_butanediol_nametrans(tmp_path):
    """Butanediol workflow — name_translation dictionary."""
    baseline = _setup(tmp_path, 'test4_butanediol_nametrans')
    name_translation = {
        'C1': 'CA', 'C2': 'CB', 'O1': 'OH',
        'H1': 'HO', 'H2': 'HA', 'H3': 'HB',
    }

    off = afm.ReadOFF(off_loc=str(SAMPLES / 'butanediol_intra.off'))
    nonbonded_tabpot = off.gen_nonbonded_tabpot()
    off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot,
                               tabpot_prefix='table',
                               name_translation=name_translation)
    bonded_tabpot = off.gen_bonded_tabpot()
    off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot, tabpot_prefix='table')
    off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top',
                               name_translation=name_translation)
    off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top',
                            name_translation=name_translation,
                            bonded_tabpot=bonded_tabpot)

    _compare_all(tmp_path, baseline)


def test5_methane_scsigma(tmp_path):
    """Methane workflow — sc_sigma for free energy calculations."""
    baseline = _setup(tmp_path, 'test5_methane_scsigma')
    sc_sigma_value = 0.3

    off = afm.ReadOFF(off_loc=str(SAMPLES / 'methane_intra.off'))
    nonbonded_tabpot = off.gen_nonbonded_tabpot(sc_sigma=sc_sigma_value)
    off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, tabpot_prefix='table')
    off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top',
                               sc_sigma=sc_sigma_value)
    off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')

    _compare_all(tmp_path, baseline)


def test6_ethane_exclpairs(tmp_path):
    """Ethane workflow — excl_pairs to exclude specific atom pairs."""
    baseline = _setup(tmp_path, 'test6_ethane_exclpairs')
    excl_pairs = [['C', 'C'], ['H', 'H']]

    off = afm.ReadOFF(off_loc=str(SAMPLES / 'ethane_intra.off'))
    nonbonded_tabpot = off.gen_nonbonded_tabpot(excl_pairs=excl_pairs)
    off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, tabpot_prefix='table')
    bonded_tabpot = off.gen_bonded_tabpot()
    off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot, tabpot_prefix='table')
    off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top',
                               excl_pairs=excl_pairs)
    off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top',
                            bonded_tabpot=bonded_tabpot)

    _compare_all(tmp_path, baseline)


def test7_methane_config(tmp_path):
    """Methane workflow — config system, charge file loading, parameter override."""
    baseline = _setup(tmp_path, 'test7_methane_config', extra_inputs=['charges.txt'])

    off = afm.ReadOFF(off_loc=str(SAMPLES / 'methane_intra.off'))
    off.set_config(
        spacing_nonbonded=0.001,
        length_nonbonded=2.5,
        tabpot_prefix='config_table',
        tabpot_dir='config_tabpot',
        scale_C6=False,
    )
    off.load_charges_from_file('charges.txt')

    nonbonded_tabpot = off.gen_nonbonded_tabpot()
    # Write using config defaults (config_tabpot/config_table_*.xvg)
    off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)
    # Write again with explicit override (override_tabpot/override_table_*.xvg)
    off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot,
                               tabpot_prefix='override_table',
                               tabpot_dir='override_tabpot')

    off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
    off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')

    _compare_all(tmp_path, baseline)


def test8_ethane_clean_workflow(tmp_path):
    """Ethane workflow — config-based API (set_config + auto-store pattern)."""
    baseline = _setup(tmp_path, 'test8_ethane_clean_workflow', extra_inputs=['charges.txt'])

    off = afm.ReadOFF(off_loc=str(SAMPLES / 'ethane_intra.off'))
    off.set_config(
        spacing_nonbonded=0.0005,
        length_nonbonded=3.0,
        spacing_bonded=0.0001,
        length_bonded=0.3,
        tabpot_prefix='table',
        tabpot_dir='',
        scale_C6=True,
    )
    off.load_charges_from_file('charges.txt')
    off.gen_nonbonded_tabpot()
    off.write_nonbonded_tabpot()
    off.gen_bonded_tabpot()
    off.write_bonded_tabpot()
    off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
    off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')

    _compare_all(tmp_path, baseline)
