"""Helper utilities for membrane-based analyses.

This module provides :class:`MembraneAnalysisBase`, an extension of
:class:`analysis.Base.AnalysisBase` with convenience methods for selecting
membrane atoms and computing density related quantities.  It encapsulates
common functionality required by multiple membrane specific analyses.
"""

import logging
import numpy as np
from .Base import AnalysisBase
import os


logger = logging.getLogger(__name__)


class MembraneAnalysisBase(AnalysisBase):
    """Common functionality for analyses that operate on membranes.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe containing the membrane system.
    lipids : list[str]
        Residue names identifying membrane phospholipids.
    NL : str
        Neutral lipid residue name (e.g., ``'TRIO'``).
    water : str
        Solvent residue name.
    start, stop, step : int, optional
        Frame range and stride passed to :class:`AnalysisBase`.
    """
    def __init__(self, universe, lipids, NL, water, start=0, stop=None, step=1):
        super().__init__(universe, start=start, stop=stop, step=step)
        self.lipids = lipids
        self.NL = NL
        self.water = water
        self.halfz = self.u.dimensions[2] / 2

    def select_near_protein(self, lipid, distance, leaflet='upper', inner_distance=None):
        """Return a selection string for lipids near the protein surface.

        Parameters
        ----------
        lipid : str
            Residue name of the lipid of interest.
        distance : float
            Cutoff distance (in Å) from the protein to include lipids.
        leaflet : {'upper', 'lower'}, optional
            Restrict selection to a single leaflet by comparing the ``z``
            coordinate to half the box height.
        inner_distance : float, optional
            If given, exclude lipids within this smaller distance from the
            protein.

        Returns
        -------
        str
            An MDAnalysis selection string.
        """
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
        """Return a list of headgroup atom names.

        Parameters
        ----------
        headgroup_atoms : list[str], optional
            Custom list of atom names to use. If ``None`` the phosphorous atom
            name ``'P'`` is returned.
        """
        return headgroup_atoms if headgroup_atoms else ["P"]

    def get_lipid_selection(self, extra_lipids=None):
        """Return an atom selection string for all membrane lipids.

        Parameters
        ----------
        extra_lipids : list[str], optional
            Additional residue names to include in the selection.

        Returns
        -------
        str
            An MDAnalysis selection that matches all requested lipids.
        """
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
        """Build an atom selection string for a single leaflet.

        Parameters
        ----------
        lipid_sel : str
            Base selection of lipid residues.
        headgroup_atoms : list[str]
            Atom names defining the headgroup region.
        tail_atoms : list[str]
            Atom names defining the tail region.
        prop : str
            Atom property used to distinguish the upper and lower leaflet,
            e.g. ``'prop z'``.
        alt : bool, optional
            If ``True`` a single marker atom defined by ``alt_marker_atom`` is
            used instead of ``headgroup_atoms``/``tail_atoms``.
        side : {'upper', 'lower'}, optional
            Select which leaflet to build the expression for.
        alt_marker_atom : str, optional
            Marker atom name to use when ``alt`` is ``True``.

        Returns
        -------
        str
            The MDAnalysis selection string.
        """

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
        """Build and store commonly used atom groups for analyses.

        Parameters
        ----------
        extra_lipids : list[str], optional
            Additional lipid residue names to include alongside ``self.lipids``.
        tail_atoms : list[str], optional
            Atom names that define the hydrophobic tails.
        headgroup_atoms : list[str], optional
            Atom names that represent the headgroup region.
        leaflet_property : str, optional
            Atom property used to determine leaflet membership.
        use_ls2, use_us2 : bool, optional
            Toggle alternative leaflet selection using a single marker atom.

        Returns
        -------
        dict
            Mapping of group labels (``'memb'``, ``'umemb'``, ``'lmemb'`` etc.)
            to :class:`MDAnalysis.core.groups.AtomGroup` objects.
        """
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
        """Determine "strong" residues based on oxygen atom counts.

        Residues are marked as strong if they contain between ``min_oxygens``
        and ``max_oxygens`` atoms whose names begin with ``strong_atom_prefix``
        and lie above the upper leaflet centre ``utz`` or below the lower
        leaflet centre ``ltz``.

        Parameters
        ----------
        trio_pos : ndarray, shape (N, 3)
            Cartesian coordinates of all neutral lipid atoms.
        utz, ltz : float
            ``z`` coordinates of the upper and lower leaflet centres.
        names : ndarray of str
            Atom names corresponding to ``trio_pos``.
        resids : ndarray of int
            Residue identifiers corresponding to ``trio_pos``.
        min_oxygens, max_oxygens : int, optional
            Bounds on the number of oxygen atoms required for a residue to be
            classified as strong.
        strong_atom_prefix : str, optional
            Prefix used to identify oxygen atoms.

        Returns
        -------
        tuple(ndarray, ndarray)
            Arrays of residue ids and counts satisfying the criteria.
        """
        if trio_pos.shape[0] != names.shape[0]:
            raise ValueError("Mismatch between trio_pos and names lengths.")

        strong_mask = ((trio_pos[:, 2] > utz) | (trio_pos[:, 2] < ltz)) & np.char.startswith(names, strong_atom_prefix)
        strong_resids, counts = np.unique(resids[strong_mask], return_counts=True)
        valid = (counts >= min_oxygens) & (counts <= max_oxygens)
        return strong_resids[valid], counts[valid]

    def density_frame(self, pos, mass, pbc, bins):
        """Return a mass density profile along ``z``.

        Parameters
        ----------
        pos : ndarray
            ``z`` coordinates of atoms.
        mass : ndarray
            Atomic masses corresponding to ``pos``.
        pbc : array-like
            Simulation box dimensions.
        bins : ndarray
            Bin edges for the histogram.

        Returns
        -------
        ndarray
            Normalised mass density per bin (g/cm^3).
        """
        dz = bins[1] - bins[0]
        hist, _ = np.histogram(pos, weights=mass, bins=bins)
        return hist / (pbc[0] * pbc[1] * dz * 0.602214)

    def calculate_overlap_and_inter(self, d0, d1, dz, threshold=0.1):
        """Compute normalized density overlap and scalar interdigitation.

        Parameters
        ----------
        d0, d1 : ndarray
            Density profiles to compare.
        dz : float
            Bin width along the ``z`` axis in ångström.
        threshold : float, optional
            Values of ``d0 + d1`` below this threshold are ignored when
            computing the overlap.

        Returns
        -------
        tuple(ndarray, float)
            Array of per-bin overlap values and the scalar interdigitation
            (integral of overlap) in nanometres.
        """
        d_sum = d0 + d1
        d_mul = d0 * d1
        mask = d_sum < threshold
        d_sum[mask] = 1
        d_mul[mask] = 0
        overlap = 4 * d_mul / d_sum**2
        interdig = np.sum(overlap) * dz / 10  # scale to nm
        return overlap, interdig

