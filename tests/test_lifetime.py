import numpy as np
from analysis.Lifetime import LifetimeAnalysis

class DummyUniverse:
    def __init__(self):
        self.dimensions = [100, 100, 40]
        self.trajectory = [None]

class DummyResidue:
    def __init__(self, resid):
        self.resid = resid


def test_finalize_calculates_lifetimes():
    u = DummyUniverse()
    analysis = LifetimeAnalysis(u, lipids=[], NL="TRIO", water="TIP3")
    analysis.trio_residues = [DummyResidue(1)]
    analysis.state = {1: [1, 1, 0, 1, 1, 1]}
    analysis._finalize()
    assert analysis.results[1] == [2, 3]
