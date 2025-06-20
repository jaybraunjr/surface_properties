from collections import defaultdict
import numpy as np
from .MembraneBase import MembraneAnalysisBase

class LifetimeAnalysis(MembraneAnalysisBase):
    def __init__(self, universe, lipids, NL, water,
                 buffer_length=20, min_oxygens=3, max_oxygens=6,
                 strong_atom_prefix="O", **kwargs):
        super().__init__(universe, lipids, NL, water, **kwargs)
        self.buffer_length = buffer_length
        self.min_oxygens = min_oxygens
        self.max_oxygens = max_oxygens
        self.strong_atom_prefix = strong_atom_prefix

    def _prepare(self):
        self.groups = self.setup_atom_groups()
        self.trio_residues = self.groups['trio'].residues
        self.state = {int(res.resid): [] for res in self.trio_residues}
        self.buffer = {int(res.resid): 0 for res in self.trio_residues}

    def _analyze_frame(self, ts):
        utz = np.mean(self.groups['umemb'].positions[:, 2])
        ltz = np.mean(self.groups['lmemb'].positions[:, 2])
        trio_pos = self.groups['trio'].positions
        names = self.groups['trio'].names.astype(str)
        resids = self.groups['trio'].resids

        strong_resids, counts = self.calculate_strong_resids(
            trio_pos, utz, ltz, names, resids,
            min_oxygens=self.min_oxygens,
            max_oxygens=self.max_oxygens,
            strong_atom_prefix=self.strong_atom_prefix
        )
        present = dict(zip(strong_resids, counts))

        for res in self.trio_residues:
            resid = int(res.resid)
            if present.get(resid, 0) >= self.min_oxygens:
                self.state[resid].append(1)
                self.buffer[resid] = self.buffer_length
            elif self.buffer[resid] > 0:
                self.state[resid].append(1)
                self.buffer[resid] -= 1
            else:
                self.state[resid].append(0)

    def _finalize(self):
        lifetimes = defaultdict(list)
        for resid, series in self.state.items():
            arr = np.array(series)
            starts = np.where(np.diff(np.pad(arr, (1, 0))) == 1)[0]
            ends = np.where(np.diff(np.pad(arr, (0, 1))) == -1)[0]
            for s, e in zip(starts, ends):
                length = e - s + 1
                if length > 1:
                    lifetimes[resid].append(length)
        self.results = dict(lifetimes)
