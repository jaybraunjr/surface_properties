import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from surface_properties.analysis.Vector import VectorOrientation

class DummyAtoms:
    def __init__(self, positions, residues=None):
        self.positions = np.array(positions)
        self.residues = residues or []

    def __len__(self):
        return len(self.positions)

    def select_atoms(self, sel):
        if '>' in sel:
            cutoff = float(sel.split('>')[1])
            pos = [p for p in self.positions if p[2] > cutoff]
        elif '<' in sel:
            cutoff = float(sel.split('<')[1])
            pos = [p for p in self.positions if p[2] < cutoff]
        else:
            pos = self.positions
        return DummyAtoms(np.array(pos))

    def center_of_mass(self):
        return np.mean(self.positions, axis=0)

def create_mock_universe(selection_map=None):
    universe = MagicMock()
    universe.trajectory = [MagicMock(frame=0)]

    selection_map = selection_map or {}

    def _select(sel):
        obj = selection_map.get(sel)
        if callable(obj) and not isinstance(obj, MagicMock):
            return obj(sel)
        return obj if obj is not None else MagicMock()

    universe.select_atoms.side_effect = _select
    return universe

def test_get_leaflet_z_cutoff_bottom():
    positions = [[0,0,-1],[0,0,0],[0,0,1],[0,0,2]]
    pl = DummyAtoms(positions)
    u = create_mock_universe({"resname POPC DOPE and name C210": pl})

    vec = VectorOrientation(u, leaflet="bottom")
    cutoff = vec._get_leaflet_z_cutoff()
    expected = np.mean([p[2] for p in positions if p[2] > np.mean([p[2] for p in positions])])
    assert cutoff == pytest.approx(expected)


def test_get_leaflet_z_cutoff_top():
    positions = [[0,0,-1],[0,0,0],[0,0,1],[0,0,2]]
    pl = DummyAtoms(positions)
    u = create_mock_universe({"resname POPC DOPE and name C210": pl})

    vec = VectorOrientation(u, leaflet="top")
    cutoff = vec._get_leaflet_z_cutoff()
    expected = np.mean([p[2] for p in positions if p[2] < np.mean([p[2] for p in positions])])
    assert cutoff == pytest.approx(expected)


def test_analyze_frame_populates_results():
    residue = MagicMock()
    residue.resid = 1

    def res_select(sel):
        if sel == "name O*":
            hg = MagicMock()
            hg.__len__.return_value = 6
            hg.center_of_mass.return_value = np.array([0,0,6])
            return hg
        else:
            tail = MagicMock()
            tail.__len__.return_value = 1
            tail.positions = np.array([[0,0,0]])
            return tail

    residue.atoms = MagicMock()
    residue.atoms.select_atoms.side_effect = res_select

    residues_group = MagicMock()
    residues_group.residues = [residue]

    u = create_mock_universe({"resname TRIO": residues_group})

    vec = VectorOrientation(u, leaflet="bottom")
    vec._prepare()

    with patch.object(vec, "_get_leaflet_z_cutoff", return_value=5.0):
        ts = MagicMock(frame=0)
        vec._analyze_frame(ts)

    assert len(vec.results["angles"]) > 0
    for tail in vec.tail_names:
        assert len(vec.results["avg_order"][tail]) == 1
        assert len(vec.results["std_order"][tail]) == 1
    assert vec.results["time_series"] == [0]


def test_analyze_frame_no_cutoff():
    residue = MagicMock()
    residue.resid = 1
    residue.atoms.select_atoms.return_value = MagicMock()
    residues_group = MagicMock()
    residues_group.residues = [residue]
    u = create_mock_universe({"resname TRIO": residues_group})

    vec = VectorOrientation(u, leaflet="bottom")
    vec._prepare()

    with patch.object(vec, "_get_leaflet_z_cutoff", return_value=None):
        vec._analyze_frame(MagicMock(frame=0))

    assert vec.results["angles"] == []
    assert vec.results["time_series"] == []
    for tail in vec.tail_names:
        assert vec.results["avg_order"][tail] == []
        assert vec.results["std_order"][tail] == []


def test_unpack_order():
    u = create_mock_universe()
    vec = VectorOrientation(u)
    vec._prepare()
    result = vec.unpack()
    assert result == (
        vec.results["angles"],
        vec.results["time_series"],
        vec.results["avg_order"],
        vec.results["std_order"],
    )
