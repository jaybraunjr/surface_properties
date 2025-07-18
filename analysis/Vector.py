"""Orientation of lipid tails relative to the membrane.

The :class:`VectorOrientation` class measures the alignment of specified tail
atoms with respect to the membrane normal.  For each frame the cosine of the
angle between the tail vector and the membrane normal is recorded along with
per-residue order parameters.
"""

from .Base import AnalysisBase
import numpy as np

class VectorOrientation(AnalysisBase):
    """Compute orientation of lipid tails with respect to the membrane normal.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Trajectory to analyse.
    residue_sel : str, optional
        Selection string for the residues of interest.
    tail_names : list[str], optional
        Names of tail atoms forming the vector pointing towards the headgroup.
    headgroup_sel : str, optional
        Selection used to locate the headgroup atoms.
    pl_selection : str, optional
        Selection of phospholipid atoms used to determine leaflet separation.
    leaflet : {'bottom', 'top'}, optional
        Which leaflet defines the interior of the membrane for headgroup
        filtering.
    expected_headgroup_count : int, optional
        Number of atoms expected for a complete headgroup selection.
    """
    def __init__(
        self,
        universe,
        residue_sel="resname TRIO",
        tail_names=None,
        headgroup_sel="name O*",
        pl_selection="resname POPC DOPE and name C210",
        leaflet="bottom",
        expected_headgroup_count=6,
        **kwargs
    ):
        super().__init__(universe, **kwargs)

        self.residue_sel = residue_sel
        self.tail_names = tail_names or ['C118', 'C218', 'C318']
        self.headgroup_sel = headgroup_sel
        self.pl_selection = pl_selection
        self.leaflet = leaflet
        self.expected_headgroup_count = expected_headgroup_count

    def _prepare(self):
        """Initialise result containers before analysis begins."""

        self.results['angles'] = []
        self.results['time_series'] = []
        self.results['avg_order'] = {tail: [] for tail in self.tail_names}
        self.results['std_order'] = {tail: [] for tail in self.tail_names}

    def _get_leaflet_z_cutoff(self):
        """Return the average z position defining the chosen leaflet."""

        pl_atoms = self.u.select_atoms(self.pl_selection)
        if len(pl_atoms) == 0:
            return None
        z_coords = pl_atoms.positions[:, 2]
        mid_z = np.mean(z_coords)
        if self.leaflet == 'bottom':
            leaflet_atoms = pl_atoms.select_atoms(f'prop z > {mid_z}')
        elif self.leaflet == 'top':
            leaflet_atoms = pl_atoms.select_atoms(f'prop z < {mid_z}')
        else:
            raise ValueError("leaflet must be 'top' or 'bottom'")
        if len(leaflet_atoms) == 0:
            return None
        return np.mean(leaflet_atoms.positions[:, 2])

    def _analyze_frame(self, ts):
        """Compute orientation metrics for one frame."""
        avg_leaflet_z = self._get_leaflet_z_cutoff()
        if avg_leaflet_z is None:
            return

        residues = self.u.select_atoms(self.residue_sel).residues
        frame_cosines = {tail: [] for tail in self.tail_names}

        for res in residues:
            headgroup = res.atoms.select_atoms(self.headgroup_sel)
            if len(headgroup) != self.expected_headgroup_count:
                continue
            com_head = headgroup.center_of_mass()
            if (self.leaflet == 'bottom' and com_head[2] <= avg_leaflet_z) or \
               (self.leaflet == 'top' and com_head[2] >= avg_leaflet_z):
                continue

            for tail in self.tail_names:
                tail_atom = res.atoms.select_atoms(f"name {tail}")
                if len(tail_atom) != 1:
                    continue
                tail_pos = tail_atom.positions[0]
                vec = com_head - tail_pos
                norm = np.linalg.norm(vec)
                if norm == 0:
                    continue
                cos_theta = np.dot(vec, [0, 0, 1]) / norm

                self.results['angles'].append({
                    'frame': ts.frame,
                    'resid': res.resid,
                    'tail_name': tail,
                    'cosine_alignment': cos_theta
                })
                frame_cosines[tail].append(cos_theta)

        for tail in self.tail_names:
            vals = frame_cosines[tail]
            self.results['avg_order'][tail].append(np.mean(vals) if vals else np.nan)
            self.results['std_order'][tail].append(np.std(vals) if vals else np.nan)

        self.results['time_series'].append(ts.frame)


    def unpack(self):
        """Return collected results as a tuple."""

        return (
            self.results.get('angles'),
            self.results.get('time_series'),
            self.results.get('avg_order'),
            self.results.get('std_order')
        )
