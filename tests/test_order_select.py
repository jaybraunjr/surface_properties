import pytest
from analysis.order import OrderParameters


class DummyUniverse:
    def __init__(self):
        self.last_selection = None
        self.trajectory = [None]

    def select_atoms(self, selection):
        self.last_selection = selection
        return "atoms"


def test_select_atoms_uses_strong_residues():
    u = DummyUniverse()
    def get_strong(ts):
        return [5, 6]
    op = OrderParameters(u, atomlists=[["C1", "H1"]], get_strong_residues=get_strong, verbose=False)
    result = op._select_atoms(None)
    assert result == "atoms"
    assert "resid 5 or resid 6" in u.last_selection


def test_select_atoms_empty_returns_none():
    u = DummyUniverse()
    def get_strong(ts):
        return []
    op = OrderParameters(u, atomlists=[["C1", "H1"]], get_strong_residues=get_strong, verbose=False)
    result = op._select_atoms(None)
    assert result is None

