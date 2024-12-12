import pytest
import numpy as np
from surface_properties.analysis.order import OrderParameters
from unittest.mock import MagicMock

# Helper function to create a mock MDAnalysis universe
def create_mock_universe():
    universe = MagicMock()
    universe.trajectory = MagicMock()
    universe.trajectory.__len__.return_value = 10  # Simulate 10 frames
    universe.atoms = MagicMock()
    universe.atoms.indices = list(range(100))  # Mock 100 atoms
    universe.dimensions = [10, 10, 10, 90, 90, 90]  # Box dimensions
    universe.select_atoms = MagicMock()
    universe.select_atoms.return_value = MagicMock()  # Mock selections
    return universe

# Test initialization of OrderParameters
def test_order_parameters_initialization():
    u = create_mock_universe()
    atomlists = [["C1", "H1", "H2"], ["C2", "H3", "H4"]]
    op = OrderParameters(u, atomlists, selection="resname POPC")
    assert len(op.C_numbers) == 2
    assert len(op.Hs_f) == 4

# Test process_atom_lists method
def test_process_atom_lists():
    u = create_mock_universe()
    atomlists = [["C1", "H1", "H2"], ["C2", "H3", "H4"]]
    op = OrderParameters(u, atomlists, selection="resname POPC")
    C_numbers, Cs, Hs_f, repeat = op.process_atom_lists()
    assert C_numbers == [1, 2]
    assert Cs == ["C1", "C2"]
    assert Hs_f == ["H1", "H2", "H3", "H4"]
    assert repeat == [2, 2]

# Test compute_OP with no strong residues
def test_compute_op_no_strong_residues():
    u = create_mock_universe()
    atomlists = [["C1", "H1"], ["C2", "H2"]]
    op = OrderParameters(u, atomlists, selection="resname POPC", start_frame=0, end_frame=10)
    u.select_atoms.return_value = MagicMock()
    u.select_atoms.return_value.residues = []  # No residues selected
    results = op.compute_OP()
    assert results.shape == (0,)

# Test for invalid frame range
def test_invalid_frame_range():
    u = create_mock_universe()
    atomlists = [["C1", "H1"]]
    with pytest.raises(ValueError, match="end_frame must be greater than start_frame"):
        OrderParameters(u, atomlists, selection="resname POPC", start_frame=10, end_frame=5)

# Test for missing selection or dynamic residues callable
def test_missing_selection_and_callable():
    u = create_mock_universe()
    atomlists = [["C1", "H1"]]
    with pytest.raises(ValueError, match="Either a selection string or a dynamic residue selection callable must be provided."):
        OrderParameters(u, atomlists)
