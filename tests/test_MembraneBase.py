from analysis.MembraneBase import MembraneAnalysisBase
import pytest
import numpy as np

class DummyUniverse:
    def __init__(self, z=40):
        self.dimensions = [100, 100, z]
        self.trajectory = [None] * 10

@pytest.fixture
def base():
    u = DummyUniverse()
    return MembraneAnalysisBase(u, lipids=[], NL='TRIO', water='TIP3')

def test_calculate_strong_resids(base):
    resid_pos = np.array([
        [0, 0, 25],
        [0, 0, 5],
        [0, 0, 15],
        [0, 0, 27]
    ])
    names = np.array(["O21", "O22", "C23", "O24"])  # only names starting with 'O' count
    resids = np.array([1, 1, 2, 3])
    utz, ltz = 20.0, 10.0

    resids_out, counts_out = base.calculate_strong_resids(
        resid_pos, utz, ltz, names, resids, min_oxygens=1, max_oxygens=5
    )

    assert set(resids_out) == {1, 3}
    assert len(counts_out) == 2
