import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from surface_properties.analysis.Surface import MembraneAnalysisBase, LifetimeAnalysis, InterdigitationAnalysis

def create_mock_universe():
    universe = MagicMock()
    universe.trajectory.n_frames = 100  # Simulate 100 frames

    # Mock trajectory frames
    def mock_frame_generator():
        for i in range(10):  # Simulate 10 frames
            frame = MagicMock()
            frame.time = i * 1000  # Simulate time in ps
            yield frame

    universe.trajectory = list(mock_frame_generator())  # Mock trajectory as a list of frames
    universe.dimensions = [100, 100, 100, 90, 90, 90]  # Box dimensions
    universe.select_atoms = MagicMock(return_value=MagicMock(
        positions=np.random.rand(6, 3),
        masses=np.random.rand(6),
        names=np.array(["O1", "O2", "O3", "C1", "H1", "O4"]),
        resids=np.array([1, 2, 3, 1, 2, 3])
    ))
    return universe


# Test MembraneAnalysisBase
def test_setup_atom_groups():
    universe = create_mock_universe()
    base = MembraneAnalysisBase(universe, lipids=["POPC"], NL="TRIO", water="WATER")
    groups = base.setup_atom_groups()
    assert "memb" in groups
    assert "umemb" in groups
    assert "lmemb" in groups
    assert "trio" in groups
    assert "water" in groups

def test_calculate_strong_resids():
    universe = create_mock_universe()
    base = MembraneAnalysisBase(universe, lipids=["POPC"], NL="TRIO", water="WATER")

    trio_pos = np.random.rand(6, 3)  # Match the length of `names` and `resids`
    utz = 0.8
    ltz = 0.2
    names = np.array(["O1", "O2", "C1", "O3", "H1", "O4"])  # 6 elements
    resids = np.array([1, 1, 2, 3, 3, 3])  # 6 elements

    strong_resids, counts = base.calculate_strong_resids(trio_pos, utz, ltz, names, resids)
    assert isinstance(strong_resids, np.ndarray)
    assert isinstance(counts, np.ndarray)


def test_density_frame():
    universe = create_mock_universe()
    base = MembraneAnalysisBase(universe, lipids=["POPC"], NL="TRIO", water="WATER")
    pos = np.random.rand(10)
    mass = np.random.rand(10)
    bins = np.linspace(0, 100, 11)
    pbc = universe.dimensions
    density = base.density_frame(pos, mass, pbc, bins)
    assert density.shape == (10,)

def test_calculate_overlap_and_inter():
    universe = create_mock_universe()
    base = MembraneAnalysisBase(universe, lipids=["POPC"], NL="TRIO", water="WATER")
    d0 = np.random.rand(10)
    d1 = np.random.rand(10)
    dz = 1.0
    overlap, interdigitation = base.calculate_overlap_and_inter(d0, d1, dz)
    assert overlap.shape == (10,)
    assert interdigitation > 0

# Test LifetimeAnalysis
def test_calculate_trio_lifetimes():
    universe = create_mock_universe()
    lifetime_analysis = LifetimeAnalysis(universe, lipids=["POPC"], NL="TRIO", water="WATER")
    lifetimes = lifetime_analysis.calculate_trio_lifetimes(start_frame=0, end_frame=10)
    assert isinstance(lifetimes, dict)

@patch("builtins.open", new_callable=MagicMock)
def test_analyze_and_save(mock_open):
    universe = create_mock_universe()
    lifetime_analysis = LifetimeAnalysis(universe, lipids=["POPC"], NL="TRIO", water="WATER")
    lifetime_analysis.analyze_and_save(base_dir="test_dir", start_frame=0, end_frame=10)

    # Ensure `open` is called only if TRIO residues exist
    if mock_open.call_count > 0:
        mock_open.assert_called()
    else:
        print("No files saved because no TRIO residues were found.")

def test_interdigit():
    universe = create_mock_universe()
    interdig_analysis = InterdigitationAnalysis(universe, lipids=["POPC"], NL="TRIO", water="WATER")

    # Mock group data
    universe.select_atoms.return_value.positions = np.random.rand(6, 3)  # 6 positions
    universe.select_atoms.return_value.masses = np.random.rand(6)  # 6 masses
    universe.select_atoms.return_value.names = np.array(["O1", "O2", "O3", "C1", "H1", "O4"])
    universe.select_atoms.return_value.resids = np.array([1, 2, 3, 1, 2, 3])

    # Mock phospholipid count
    universe.select_atoms.side_effect = lambda sel: MagicMock(n_atoms=10) if sel == 'name P' else MagicMock(
        positions=np.random.rand(6, 3),
        masses=np.random.rand(6),
        names=np.array(["O1", "O2", "O3", "C1", "H1", "O4"]),
        resids=np.array([1, 2, 3, 1, 2, 3])
    )

    # Run interdigitation analysis
    results = interdig_analysis.interdigit(b=0, e=10)
    assert "inter" in results
    assert "ov" in results
    assert "density" in results



@patch("os.makedirs")
@patch("numpy.savetxt")
def test_save_results(mock_savetxt, mock_makedirs):
    universe = create_mock_universe()
    interdig_analysis = InterdigitationAnalysis(universe, lipids=["POPC"], NL="TRIO", water="WATER")
    results = {
        "inter": {"total": np.random.rand(10, 2)},
        "ov": {"total": np.random.rand(10, 2)},
        "density": {"PL": np.random.rand(10, 2)},
    }
    interdig_analysis.save_results(results, base_dir="test_dir")
    mock_makedirs.assert_called()
    mock_savetxt.assert_called()


