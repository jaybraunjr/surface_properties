import pytest
from analysis.Base import AnalysisBase

class DummyUniverse:
    def __init__(self, n_frames=10):
        self.trajectory = [object() for _ in range(n_frames)]

class RecordingAnalysis(AnalysisBase):
    def _prepare(self):
        self.calls = 0

    def _analyze_frame(self, ts):
        self.calls += 1


def test_run_slice():
    u = DummyUniverse()
    analysis = RecordingAnalysis(u, start=2, stop=8, step=2, verbose=False)
    analysis.run()
    assert analysis.calls == 3
    assert analysis._frame_index == 3
