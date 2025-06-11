from .base import AnalysisBase  # your custom base
import numpy as np
import logging
from unittest.mock import MagicMock

class OrderParameters(AnalysisBase):
    def __init__(self, u, atomlists, selection=None, get_strong_residues=None,
                 start_frame=0, end_frame=None, **kwargs):
        if selection is None and get_strong_residues is None:
            raise ValueError("Either a selection string or a dynamic residue selection callable must be provided.")
        super().__init__(u, start=start_frame, stop=end_frame, **kwargs)
        self.atomlists = atomlists
        self.selection = selection
        self.get_strong_residues = get_strong_residues
        self.C_numbers, self.Cs, self.Hs_f, self.repeat = self.process_atom_lists()

    def process_atom_lists(self):
        C_numbers, Cs, Hs, repeat = [], [], [], []
        for atoms in self.atomlists:
            # C_number = ''.join(filter(str.isdigit, atoms[0]))
            # C_numbers.append(int(C_number))
            C_numbers.append(len(C_numbers) + 1)  # sequential: 1, 2, 3, ...
            Cs.append(atoms[0])
            Hs.extend(atoms[1:])
            repeat.append(len(atoms) - 1)
        return C_numbers, Cs, Hs, repeat

    def _prepare(self):
        self.results['order_parameters'] = []

    def _select_atoms(self, ts):
        if self.get_strong_residues:
            strong_residues = self.get_strong_residues(ts)
            if not strong_residues:
                return None
            sel = "resid " + " or resid ".join(map(str, strong_residues))
            return self.u.select_atoms(sel)
        elif self.selection:
            return self.u.select_atoms(self.selection)
        return None

    def _analyze_frame(self, ts):
        atoms = self._select_atoms(ts)
        if atoms is None or len(atoms) == 0:
            self.results['order_parameters'].append(np.zeros(len(self.Cs)))
            return

        valid_C, valid_H = [], []
        for res in atoms.residues:
            Cs = res.atoms.select_atoms("name " + " ".join(self.Cs))
            Hs = res.atoms.select_atoms("name " + " ".join(self.Hs_f))
            if len(Cs) and len(Hs):
                valid_C.extend(Cs.indices)
                valid_H.extend(Hs.indices)

        if not valid_C or not valid_H:
            self.results['order_parameters'].append(np.zeros(len(self.Cs)))
            return

        group1 = self.u.atoms[valid_C]
        group2 = self.u.atoms[valid_H]
        natoms = len(self.Cs)
        nmols = len(group1.positions) // natoms
        repeats = self.repeat * nmols

        p1 = np.repeat(group1.positions, repeats, axis=0)
        p2 = group2.positions
        dp = p2 - p1
        norm = np.linalg.norm(dp, axis=-1)
        norm[norm == 0] = 1e-8
        cos_theta = dp[:, 2] / norm
        S = -0.5 * (3 * cos_theta**2 - 1)

        new_S = self._average_over_hydrogens(S, repeats)
        new_S.shape = (nmols, natoms)
        self.results['order_parameters'].append(np.mean(new_S, axis=0))

    def _finalize(self):
        arr = np.array(self.results['order_parameters'])
        self.results['average'] = np.mean(arr, axis=0)
        self.results['output'] = np.transpose([self.C_numbers, self.results['average']])

    def _average_over_hydrogens(self, x, reps):
        out, i = [], 0
        for rep in reps:
            out.append(np.mean(x[i:i+rep]))
            i += rep
        return np.array(out)

    def compute_OP(self):
        self._prepare()
        end = self.stop or len(self.u.trajectory)
        for frame in range(self.start, end, self.step):
            self._analyze_frame(MagicMock(frame=frame))
        self._finalize()
        return np.array(self.results.get('order_parameters'))

