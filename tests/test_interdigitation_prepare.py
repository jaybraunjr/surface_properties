import numpy as np
import pytest
from analysis.Surface import InterdigitationAnalysis

class DummyGroup:
    def __init__(self, n_atoms=0):
        self.n_atoms = n_atoms
        self.positions = np.zeros((n_atoms, 3))
        self.masses = np.ones(n_atoms)
        self.names = np.array([])
        self.resids = np.array([])

class DummyUniverse:
    def __init__(self, z=40):
        self.dimensions = [100, 100, z]
        self.trajectory = [None]
    def select_atoms(self, selection):
        return DummyGroup()

def test_interdigitation_prepare():
    z = 30
    u = DummyUniverse(z=z)
    analysis = InterdigitationAnalysis(u, lipids=[], NL='TRIO', water='TIP3')
    analysis._prepare()

    assert len(analysis.bins) == analysis.nbins + 1
    assert pytest.approx(analysis.dz) == z / analysis.nbins

    for key in ["memb", "trio", "water", "umemb", "lmemb"]:
        assert key in analysis.groups
