import numpy as np
from .Base import AnalysisBase
from .MembraneBase import MembraneAnalysisBase
import os


# Interdigitation analysis class. Inherits from AnalysisBase in Base.py.
# This impliments the membrane logic.


class InterdigitationAnalysis(AnalysisBase):
    
    """
    Analysis class for computing interdigitation and density overlap
    between membrane lipids and neutral lipids (NL) such as TRIO.
    
    Runs per-frame analysis and aggregates results, supporting output
    of strong residue lists and summary statistics.
    """
    def __init__(self, universe, lipids, NL, water,
                 start=0, stop=None, step=1,
                 strong_resid_list_name='strong_resid_list',
                 extra_lipids=None, tail_atoms=None, headgroup_atoms=None,
                 leaflet_property="prop z", use_ls2=False, use_us2=False,
                 min_oxygens=3, max_oxygens=6, strong_atom_prefix="O"):
        
        """
        Initialize the InterdigitationAnalysis.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            The MDAnalysis Universe object for the simulation.
        lipids : list of str
            List of membrane lipid residue names.
        NL : str
            Neutral lipid residue name (e.g., 'TRIO').
        water : str
            Water residue name (e.g., 'WAT').
        start, stop, step : int, optional
            Frame range and stride for analysis.
        strong_resid_list_name : str, optional
            Base name for saving strong residue lists.
        extra_lipids, tail_atoms, headgroup_atoms : list, optional
            Additional selection options.
        leaflet_property : str, optional
            Atom property to use for leaflet selection.
        use_ls2, use_us2 : bool, optional
            Alternative leaflet selection logic.
        min_oxygens, max_oxygens : int, optional
            Criteria for strong residue identification.
        strong_atom_prefix : str, optional
            Atom name prefix (default "O") for strong residue selection.
        """

        super().__init__(universe, start=start, stop=stop, step=step)

        self.base = MembraneAnalysisBase(universe, lipids, NL, water)
        self.strong_resid_list_name = strong_resid_list_name
        self.extra_lipids = extra_lipids
        self.tail_atoms = tail_atoms
        self.headgroup_atoms = headgroup_atoms
        self.leaflet_property = leaflet_property
        self.use_ls2 = use_ls2
        self.use_us2 = use_us2
        self.min_oxygens = min_oxygens
        self.max_oxygens = max_oxygens
        self.strong_atom_prefix = strong_atom_prefix
        self.groups = None
        self.nbins = 50
        self.bins = None
        self.dz = None
        self.numP = None

        self.frame_data = []

    def _prepare(self):
        self.tail_atoms = self.base.get_tail_atoms(self.tail_atoms)
        self.groups = self.base.setup_atom_groups(
            extra_lipids=self.extra_lipids,
            tail_atoms=self.tail_atoms,
            headgroup_atoms=self.headgroup_atoms,
            leaflet_property=self.leaflet_property,
            use_ls2=self.use_ls2,
            use_us2=self.use_us2
        )
        self.numP = self.u.select_atoms("name P").n_atoms
        pbc = self.u.dimensions
        self.bins = np.linspace(0, pbc[2], self.nbins + 1)
        self.dz = self.bins[1] - self.bins[0]

    def _analyze_frame(self, ts):
        utz = np.mean(self.groups['umemb'].positions[:, 2])
        ltz = np.mean(self.groups['lmemb'].positions[:, 2])

        nl_key = self.base.NL.lower()
        nl_group = self.groups[nl_key]
        nl_pos = nl_group.positions
        names = nl_group.names.astype(str)
        resids = nl_group.resids

        strong_resids, _ = self.base.calculate_strong_resids(
            nl_pos, utz, ltz, names, resids,
            min_oxygens=self.min_oxygens, max_oxygens=self.max_oxygens,
            strong_atom_prefix=self.strong_atom_prefix
        )

        d0 = self.base.density_frame(
            self.groups['memb'].select_atoms(f"name {' '.join(self.tail_atoms)}").positions[:, 2],
            self.groups['memb'].select_atoms(f"name {' '.join(self.tail_atoms)}").masses,
            self.u.dimensions, self.bins
        )
        d1 = self.base.density_frame(
            nl_group.positions[:, 2],
            nl_group.masses,
            self.u.dimensions, self.bins
        )
        d2 = self.base.density_frame(
            nl_group.positions[np.isin(resids, strong_resids), 2],
            nl_group.masses[np.isin(resids, strong_resids)],
            self.u.dimensions, self.bins
        )
        d3 = self.base.density_frame(
            nl_group.positions[np.isin(resids, strong_resids, invert=True), 2],
            nl_group.masses[np.isin(resids, strong_resids, invert=True)],
            self.u.dimensions, self.bins
        )
        d_water = self.base.density_frame(
            self.groups['water'].positions[:, 2],
            self.groups['water'].masses,
            self.u.dimensions, self.bins
        )

        tov, tin = self.base.calculate_overlap_and_inter(d0, d1, self.dz)
        sov, sin = self.base.calculate_overlap_and_inter(d0, d2, self.dz)
        wov, win = self.base.calculate_overlap_and_inter(d0, d3, self.dz)

        nl_to_pl = len(strong_resids) / self.numP if self.numP > 0 else 0

        self.frame_data.append({
            'time': ts.time / 1000,
            'z': self.u.dimensions[2] / 10,
            'strong_resids': strong_resids.tolist(),
            'inter': (tin, sin, win),
            'ov': (tov, sov, wov),
            'densities': (d0, d1, d2, d3, d_water),
            'num_strong': len(strong_resids),
            'nl_to_pl': nl_to_pl
        })

    def _finalize(self):
        times = np.array([d['time'] for d in self.frame_data])
        zs = np.array([d['z'] for d in self.frame_data])
        XX = np.linspace(0, np.mean(zs), self.nbins)

        def average_stack(key, idx):
            return np.mean([d[key][idx] for d in self.frame_data], axis=0)
        
        NL_label = self.base.NL

        self.results = {
            'inter': {
                'total': np.transpose([times, [d['inter'][0] for d in self.frame_data]]),
                'strong': np.transpose([times, [d['inter'][1] for d in self.frame_data]]),
                'weak': np.transpose([times, [d['inter'][2] for d in self.frame_data]])
            },
            'ov': {
                'total': np.transpose([XX, average_stack('ov', 0)]),
                'strong': np.transpose([XX, average_stack('ov', 1)]),
                'weak': np.transpose([XX, average_stack('ov', 2)])
            },
            'ratio': {
                'num': np.transpose([times, [d['num_strong'] for d in self.frame_data]]),
                f'{NL_label.lower()}-to-pl': np.transpose([times, [d['nl_to_pl'] for d in self.frame_data]])
            },
            'density': {
                'PL': np.transpose([XX, average_stack('densities', 0)]),
                NL_label: np.transpose([XX, average_stack('densities', 1)]),
                'water': np.transpose([XX, average_stack('densities', 4)])
            },
            'strong_residues': [d['strong_resids'] for d in self.frame_data]
        }

        with open(f"{self.strong_resid_list_name}.txt", "w") as f:
            f.write('\n'.join(map(str, self.results['strong_residues'])))



    def save_results(self, out_dir='results', prefix=None):

        os.makedirs(out_dir, exist_ok=True)
        base_name = prefix or self.strong_resid_list_name

        strong_res_path = os.path.join(out_dir, f"{base_name}_strong_residues.txt")
        with open(strong_res_path, "w") as f:
            for frame_resids in self.results['strong_residues']:
                f.write(f"{frame_resids}\n")
        print(f"Strong residue lists saved to {strong_res_path}")

