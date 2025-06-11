import numpy as np
import logging
from .Base import AnalysisBase
from tqdm import tqdm  
from collections import defaultdict
import numpy as np
import os

import matplotlib.pyplot as plt
import MDAnalysis as mda



logger = logging.getLogger(__name__)


class MembraneAnalysisBase(AnalysisBase):
    def __init__(self, universe, lipids, NL, water, start=0, stop=None, step=1):
        super().__init__(universe, start=start, stop=stop, step=step)
        self.lipids = lipids
        self.NL = NL
        self.water = water
        self.halfz = self.u.dimensions[2] / 2


    def get_tail_atoms(self, tail_atoms=None):
        if tail_atoms is None:
            tail_atoms = [f"C2{i}" for i in range(2, 22)] + [f"C3{i}" for i in range(2, 22)]
        return tail_atoms

    def get_headgroup_atoms(self, headgroup_atoms=None):
        return headgroup_atoms if headgroup_atoms else ["P"]

    def get_lipid_selection(self, extra_lipids=None):
        selected_lipids = self.lipids + (extra_lipids or [])
        return " or ".join(f"resname {lipid}" for lipid in selected_lipids)

    def get_leaflet_selection(self, lipid_sel, headgroup_atoms, tail_atoms, prop, alt=False, side='upper'):
        z_cmp = ">" if side == 'upper' else "<"
        if alt:
            return f"(same residue as ({lipid_sel}) and name C24 and {prop}{z_cmp}{self.halfz})"
        else:
            heads = " ".join(headgroup_atoms)
            tails = " ".join(tail_atoms)
            return f"(same residue as ({lipid_sel}) and name {heads} and {prop}{z_cmp}{self.halfz}) and name {tails}"

    def setup_atom_groups(self, extra_lipids=None, tail_atoms=None, headgroup_atoms=None,
                          leaflet_property="prop z", use_ls2=False, use_us2=False):
        lipid_sel = self.get_lipid_selection(extra_lipids)
        tail_atoms = self.get_tail_atoms(tail_atoms)
        headgroup_atoms = self.get_headgroup_atoms(headgroup_atoms)

        upper_sel = self.get_leaflet_selection(lipid_sel, headgroup_atoms, tail_atoms, leaflet_property, alt=use_us2, side='upper')
        lower_sel = self.get_leaflet_selection(lipid_sel, headgroup_atoms, tail_atoms, leaflet_property, alt=use_ls2, side='lower')

        groups = {
            "memb": self.u.select_atoms(lipid_sel),
            "umemb": self.u.select_atoms(upper_sel),
            "lmemb": self.u.select_atoms(lower_sel),
            "trio": self.u.select_atoms(f"resname {self.NL}"),
            "water": self.u.select_atoms(f"resname {self.water}")
        }

        logger.debug(f"Tail atoms used: {tail_atoms}")
        logger.debug(f"Membrane selection: {lipid_sel}")
        logger.debug(f"Upper leaflet selection: {upper_sel}")
        logger.debug(f"Lower leaflet selection: {lower_sel}")

        return groups

    def calculate_strong_resids(self, trio_pos, utz, ltz, names, resids,
                                min_oxygens=3, max_oxygens=6, strong_atom_prefix="O"):
        if trio_pos.shape[0] != names.shape[0]:
            raise ValueError("Mismatch between trio_pos and names lengths.")

        strong_mask = []
        for pos, name in zip(trio_pos, names):
            flag = (pos[2] > utz or pos[2] < ltz) and name.startswith(strong_atom_prefix)
            strong_mask.append(flag)
        selected = [r for r, m in zip(resids, strong_mask) if m]
        count_map = {}
        for r in selected:
            count_map[r] = count_map.get(r, 0) + 1
        strong_resids = []
        counts = []
        for resid, count in count_map.items():
            if min_oxygens <= count <= max_oxygens:
                strong_resids.append(resid)
                counts.append(count)
        return np.array(strong_resids), np.array(counts)

    def density_frame(self, pos, mass, pbc, bins):
        dz = bins[1] - bins[0]
        hist, _ = np.histogram(pos, weights=mass, bins=bins)
        return hist / (pbc[0] * pbc[1] * dz * 0.602214)

    def calculate_overlap_and_inter(self, d0, d1, dz, threshold=0.1):
        d_sum = d0 + d1
        d_mul = d0 * d1
        mask = d_sum < threshold
        d_sum[mask] = 1
        d_mul[mask] = 0
        overlap = 4 * d_mul / d_sum**2
        interdig = np.sum(overlap) * dz / 10  # scale to nm
        return overlap, interdig
