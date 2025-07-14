import logging
import numpy as np
from .Base import AnalysisBase
import os

# A helper class for membrane analyses.


logger = logging.getLogger(__name__)


class MembraneAnalysisBase(AnalysisBase):
    def __init__(self, universe, lipids, NL, water, start=0, stop=None, step=1):
        super().__init__(universe, start=start, stop=stop, step=step)
        self.lipids = lipids
        self.NL = NL
        self.water = water
        self.halfz = self.u.dimensions[2] / 2

    def select_near_protein(self, lipid, distance, leaflet='upper', inner_distance=None):
        halfz = self.u.dimensions[2]/2
        z_filter = f"prop z > {halfz}" if leaflet == 'upper' else f"prop z < {halfz}"
        near = f"resname {lipid} and (around {distance} protein) and {z_filter}"
        if inner_distance is not None:
            near = f"{near} and not (around {inner_distance} protein)"
        return near
    

    def get_tail_atoms(self, tail_atoms=None):
        """Return a list of tail atom names (default is C2* and C3*)."""
        if tail_atoms is None:
            tail_atoms = [f"C2{i}" for i in range(2, 22)] + [f"C3{i}" for i in range(2, 22)]
        return tail_atoms

    def get_headgroup_atoms(self, headgroup_atoms=None):
        return headgroup_atoms if headgroup_atoms else ["P"]

    def get_lipid_selection(self, extra_lipids=None):
        """Return atom selection string for all membrane lipids."""
        selected_lipids = self.lipids + (extra_lipids or [])
        return " or ".join(f"resname {lipid}" for lipid in selected_lipids)

    # def get_leaflet_selection(self, lipid_sel, headgroup_atoms, tail_atoms, prop, alt=False, side='upper'):
    #     """Build a leaflet selection string."""
    #     z_cmp = ">" if side == 'upper' else "<"
    #     if alt:
    #         return f"(same residue as ({lipid_sel}) and name C24 and {prop}{z_cmp}{self.halfz})"
    #     else:
    #         heads = " ".join(headgroup_atoms)
    #         tails = " ".join(tail_atoms)
    #         return f"(same residue as ({lipid_sel}) and name {heads} and {prop}{z_cmp}{self.halfz}) and name {tails}"
        
    def get_leaflet_selection(
        self,
        lipid_sel,
        headgroup_atoms,
        tail_atoms,
        prop,
        alt=False,
        side='upper',
        alt_marker_atom=None
    ):

        z_cmp = ">" if side == 'upper' else "<"
        if alt:
            marker = alt_marker_atom if alt_marker_atom else "C24"
            return f"(same residue as ({lipid_sel}) and name {marker} and {prop}{z_cmp}{self.halfz})"
        else:
            heads = " ".join(headgroup_atoms)
            tails = " ".join(tail_atoms)
            return (
                f"(same residue as ({lipid_sel}) and name {heads} and {prop}{z_cmp}{self.halfz}) "
                f"and name {tails}"
            )


    def setup_atom_groups(self, extra_lipids=None, tail_atoms=None, headgroup_atoms=None,
                          leaflet_property="prop z", use_ls2=False, use_us2=False):
        """Main selector for atom groups used in density and overlap calculations."""
        lipid_sel = self.get_lipid_selection(extra_lipids)
        tail_atoms = self.get_tail_atoms(tail_atoms)
        headgroup_atoms = self.get_headgroup_atoms(headgroup_atoms)

        upper_sel = self.get_leaflet_selection(lipid_sel, headgroup_atoms, tail_atoms, leaflet_property, alt=use_us2, side='upper')
        lower_sel = self.get_leaflet_selection(lipid_sel, headgroup_atoms, tail_atoms, leaflet_property, alt=use_ls2, side='lower')

        groups = {
            "memb": self.u.select_atoms(lipid_sel),
            "umemb": self.u.select_atoms(upper_sel),
            "lmemb": self.u.select_atoms(lower_sel),
            self.NL.lower(): self.u.select_atoms(f"resname {self.NL}"),
            "water": self.u.select_atoms(f"resname {self.water}")
        }

        logger.debug(f"Tail atoms used: {tail_atoms}")
        logger.debug(f"Membrane selection: {lipid_sel}")
        logger.debug(f"Upper leaflet selection: {upper_sel}")
        logger.debug(f"Lower leaflet selection: {lower_sel}")

        return groups

    def calculate_strong_resids(self, trio_pos, utz, ltz, names, resids,
                                min_oxygens=3, max_oxygens=6, strong_atom_prefix="O"):
        """Determine 'strong' residues based on oxygen atom counts above/below leaflet centers."""
        if trio_pos.shape[0] != names.shape[0]:
            raise ValueError("Mismatch between trio_pos and names lengths.")

        strong_mask = ((trio_pos[:, 2] > utz) | (trio_pos[:, 2] < ltz)) & np.char.startswith(names, strong_atom_prefix)
        strong_resids, counts = np.unique(resids[strong_mask], return_counts=True)
        valid = (counts >= min_oxygens) & (counts <= max_oxygens)
        return strong_resids[valid], counts[valid]

    def density_frame(self, pos, mass, pbc, bins):
        """Return mass density profile in z."""
        dz = bins[1] - bins[0]
        hist, _ = np.histogram(pos, weights=mass, bins=bins)
        return hist / (pbc[0] * pbc[1] * dz * 0.602214)

    def calculate_overlap_and_inter(self, d0, d1, dz, threshold=0.1):
        """Compute normalized density overlap and scalar interdigitation."""
        d_sum = d0 + d1
        d_mul = d0 * d1
        mask = d_sum < threshold
        d_sum[mask] = 1
        d_mul[mask] = 0
        overlap = 4 * d_mul / d_sum**2
        interdig = np.sum(overlap) * dz / 10  # scale to nm
        return overlap, interdig

