"""
Parsing unit tests for afmtogmx.

Each test parses a sample .off file and compares off.nonbonded and off.bonded
against stored JSON baselines in test/baseline_outputs/parsing/<stem>/.

This catches bugs in the .off parser (functions.py) independently of whether
those bugs happen to cancel out in the final .xvg/.top outputs.

Tuple keys in the data structures are serialized as their string representation
(e.g., "('H', 'H')") so they survive JSON round-trips.
"""
import sys
import json
import pathlib

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))
import afmtogmx as afm

TEST_DIR = pathlib.Path(__file__).parent
SAMPLES  = TEST_DIR / 'sample_off_files'
BASELINE = TEST_DIR / 'baseline_outputs' / 'parsing'

OFF_FILES = [
    'methane_intra.off',
    'ethane_intra.off',
    'water_intra.off',
    'butanediol_intra.off',
    'big_alanine.off',
    'curcubituril.off',
]


def _normalize(obj):
    """Recursively convert to JSON-serializable form.

    Tuple keys become their repr string; tuples in values become lists.
    """
    if isinstance(obj, dict):
        return {str(k): _normalize(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_normalize(item) for item in obj]
    else:
        return obj


def _load_baseline(stem: str, attr: str) -> dict:
    path = BASELINE / stem / f'{attr}.json'
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Parametrized parsing tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('off_file', OFF_FILES)
def test_nonbonded_parsing(off_file):
    """Parsed off.nonbonded matches stored baseline for each sample .off file."""
    stem = off_file.replace('.off', '')
    off  = afm.ReadOFF(off_loc=str(SAMPLES / off_file))

    actual   = _normalize(off.nonbonded)
    baseline = _load_baseline(stem, 'nonbonded')

    assert actual == baseline, (
        f"off.nonbonded mismatch for {off_file}\n"
        f"  actual keys:   {sorted(actual.keys())}\n"
        f"  baseline keys: {sorted(baseline.keys())}"
    )


@pytest.mark.parametrize('off_file', OFF_FILES)
def test_bonded_parsing(off_file):
    """Parsed off.bonded matches stored baseline for each sample .off file."""
    stem = off_file.replace('.off', '')
    off  = afm.ReadOFF(off_loc=str(SAMPLES / off_file))

    actual   = _normalize(off.bonded)
    baseline = _load_baseline(stem, 'bonded')

    assert actual == baseline, (
        f"off.bonded mismatch for {off_file}\n"
        f"  actual molecules:   {sorted(actual.keys())}\n"
        f"  baseline molecules: {sorted(baseline.keys())}"
    )
