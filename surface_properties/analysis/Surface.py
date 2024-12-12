from collections import defaultdict
import numpy as np
import os
import json
from tqdm import tqdm



class MembraneAnalysisBase:
    def __init__(self, universe, lipids, NL, water):
        self.u = universe
        self.lipids = lipids
        self.NL = NL
        self.water = water

    def setup_atom_groups(self, use_ls2=False, use_us2=False):
        halfz = self.u.dimensions[2] / 2
        C2 = ' '.join([f'C2{i}' for i in range(2, 22)])
        C3 = ' '.join([f'C3{i}' for i in range(2, 22)])
        tail_atoms = ' '.join(C2.split() + C3.split())

        lipid_selection = ' or '.join([f'resname {lipid}' for lipid in self.lipids])
        us = f'(same residue as ({lipid_selection}) and name P and prop z>{halfz}) and (name {tail_atoms})'
        ls = f'(same residue as ({lipid_selection}) and name P and prop z<{halfz}) and (name {tail_atoms})'

        us2 = f'(same residue as ({lipid_selection}) and name C24 and prop z>{halfz})'
        ls2 = f'(same residue as ({lipid_selection}) and name C24 and prop z<{halfz})'

        # Use instance attributes to select upper and lower leaflets
        upper_selection = us2 if use_us2 else us
        lower_selection = ls2 if use_ls2 else ls

        groups = {
            'memb': self.u.select_atoms(lipid_selection),
            'umemb': self.u.select_atoms(upper_selection),
            'lmemb': self.u.select_atoms(lower_selection),
            'trio': self.u.select_atoms(f'resname {self.NL}'),
            'water': self.u.select_atoms(f'resname {self.water}')
        }
        return groups
    
    # def calculate_strong_resids(self, trio_pos, utz, ltz, names, resids, min_oxygens=3, max_oxygens=6):
    #     boolArray = ((trio_pos[:, 2] > utz) | (trio_pos[:, 2] < ltz)) & np.char.startswith(names, 'O')
    #     strong_resids = resids[boolArray]
    #     r = np.unique(strong_resids, return_counts=True)
    #     boolArray = (r[1] >= min_oxygens) & (r[1] <= max_oxygens)
    #     strong_resids = r[0][boolArray]
    #     counts = r[1][boolArray]
    #     return strong_resids, counts
    
    def calculate_strong_resids(self, trio_pos, utz, ltz, names, resids, min_oxygens=3, max_oxygens=6):
        if len(trio_pos) != len(names):
            raise ValueError("Mismatch between trio_pos and names lengths.")
            
        boolArray = ((trio_pos[:, 2] > utz) | (trio_pos[:, 2] < ltz)) & np.char.startswith(names, 'O')
        strong_resids = resids[boolArray]
        r = np.unique(strong_resids, return_counts=True)
        boolArray = (r[1] >= min_oxygens) & (r[1] <= max_oxygens)
        strong_resids = r[0][boolArray]
        counts = r[1][boolArray]
        return strong_resids, counts

    

    def density_frame(self, pos, mass, pbc, bins):
        dz = bins[1] - bins[0]
        h, _ = np.histogram(pos, weights=mass, bins=bins)
        h /= pbc[0] * pbc[1] * dz * 0.602214
        return h

    def calculate_overlap_and_inter(self, d0, d1, dz, threshold=0.1):
        d_sum = d0 + d1
        d_mul = d0 * d1
        d_sum[d_sum < threshold] = 1
        d_mul[d_sum < threshold] = 0
        overlap = 4 * d_mul / d_sum**2
        interdigitation = np.sum(overlap) * dz / 10  
        return overlap, interdigitation



class LifetimeAnalysis(MembraneAnalysisBase):
    def __init__(self, universe, lipids, NL, water, use_ls2=False, use_us2=False, buffer_frames=1, min_oxygens=3, buffer_length=20):
        super().__init__(universe, lipids, NL, water)
        self.use_ls2 = use_ls2
        self.use_us2 = use_us2
        self.buffer_frames = buffer_frames
        self.min_oxygens = min_oxygens  # Minimum number of oxygens required
        self.buffer_length = buffer_length 

    def calculate_trio_lifetimes(self, start_frame=0, end_frame=None, step_frame=1):
        if end_frame is None:
            end_frame = self.u.trajectory.n_frames
        groups = self.setup_atom_groups(use_ls2=self.use_ls2, use_us2=self.use_us2)
        TRIO_residues = groups['trio'].residues

        if len(TRIO_residues) == 0:
            print("No TRIO residues found!")
            return {}

        TRIO_states = {int(res.resid): [] for res in TRIO_residues}
        buffer_counters = {int(res.resid): 0 for res in TRIO_residues}

        trajectory = tqdm(self.u.trajectory[start_frame:end_frame:step_frame], desc='Processing frames', unit='frame')

        for ts in trajectory:
            z_center = np.mean(groups['memb'].positions[:, 2])
            upper_lipid_atoms = groups['umemb']
            lower_lipid_atoms = groups['lmemb']
            utz = np.mean(upper_lipid_atoms.positions[:, 2])
            ltz = np.mean(lower_lipid_atoms.positions[:, 2])
            trio_pos = groups['trio'].positions
            names = groups['trio'].names.astype(str)
            resids = groups['trio'].resids
            strong_resids, counts = self.calculate_strong_resids(trio_pos, utz, ltz, names, resids)
            counts_dict = dict(zip(strong_resids, counts))

            for res in TRIO_residues:
                resid = int(res.resid)
                count = counts_dict.get(resid, 0)
                if count >= self.min_oxygens:
                    TRIO_states[resid].append(1)
                    buffer_counters[resid] = self.buffer_length  # Reset buffer counter to buffer_length
                else:
                    if buffer_counters[resid] > 0:
                        TRIO_states[resid].append(1)
                        buffer_counters[resid] -= 1
                    else:
                        TRIO_states[resid].append(0)

        lifetimes = defaultdict(list)
        for resid, state_list in TRIO_states.items():
            state_array = np.array(state_list)
            changes = np.diff(state_array)
            start_indices = np.where(changes == 1)[0] + 1
            end_indices = np.where(changes == -1)[0] + 1

            if state_array[0] == 1:
                start_indices = np.insert(start_indices, 0, 0)
            if state_array[-1] == 1:
                end_indices = np.append(end_indices, len(state_array))

            for start, end in zip(start_indices, end_indices):
                lifetime = (end - start)
                if lifetime > 1:  # Skip lifetimes that are only 1 frame
                    lifetimes[int(resid)].append(lifetime)

        return lifetimes

    def analyze_and_save(self, base_dir, start_frame=0, end_frame=None, step_frame=1):
        lifetimes = self.calculate_trio_lifetimes(start_frame=start_frame, end_frame=end_frame, step_frame=step_frame)
        if len(lifetimes) == 0:
            print("No lifetimes to save.")
            return
        os.makedirs(base_dir, exist_ok=True)
        filename = os.path.join(base_dir, 'trio_lifetimes.json')
        json_lifetimes = {str(resid): [int(lifetime) for lifetime in lifetimes[resid]] for resid in lifetimes}
        with open(filename, 'w') as f:
            json.dump(json_lifetimes, f)

        print(f"Saved TRIO lifetimes to {filename}")


class InterdigitationAnalysis(MembraneAnalysisBase):
    def __init__(self, universe, lipids, NL, water, strong_resid_list_name='strong_resid_list'):
        super().__init__(universe, lipids, NL, water)
        self.strong_resid_list_name = strong_resid_list_name

    def calculate_densities(self, groups, pbc, bins):
        # Calculate density profiles for each component
        d0 = self.density_frame(groups['memb'].positions[:, 2], groups['memb'].masses, pbc=pbc, bins=bins)
        d1 = self.density_frame(groups['trio'].positions[:, 2], groups['trio'].masses, pbc=pbc, bins=bins)
        d_water = self.density_frame(groups['water'].positions[:, 2], groups['water'].masses, pbc=pbc, bins=bins)
        return d0, d1, d_water


    def calculate_densities_for_resids(self, trio_pos, resids, strong_resids, groups, pbc, bins, invert=False):
        boolArray = np.isin(resids, strong_resids, invert=invert)
        pp = trio_pos[boolArray]
        mm = groups['trio'].masses[boolArray]
        return self.density_frame(pp[:, 2], mm, pbc, bins)


    def interdigit(self, nbins=50, nblocks=5, b=0, e=None):
        times, zs, total_inter, strong_inter, weak_inter = ([] for _ in range(5))
        d0_densities, d1_densities, d2_densities, d3_densities, d_water_series = ([] for _ in range(5))
        total_ov, strong_ov, weak_ov, strong_num = ([] for _ in range(4))
        strong_residue_ids = [] # List of strong residue IDs
        groups = self.setup_atom_groups()
        names = groups['trio'].names.astype(str)
        resids = groups['trio'].resids
        numP = self.u.select_atoms('name P').n_atoms

        # Loop over frames
        trajectory = tqdm(self.u.trajectory[b:e], desc='Processing frames', unit='frame')
        for ts in trajectory:
            pbc = self.u.dimensions
            bins = np.linspace(0, pbc[2], nbins + 1)
            dz = bins[1] - bins[0]
            trio_pos = groups['trio'].positions
            utz = np.mean(groups['umemb'].positions[:, 2])
            ltz = np.mean(groups['lmemb'].positions[:, 2])
            d0, d1, d_water = self.calculate_densities(groups, pbc, bins)
            d0_densities.append(d0)
            d1_densities.append(d1)
            d_water_series.append(d_water)
            tov, tin = self.calculate_overlap_and_inter(d0, d1, dz)
            total_ov.append(tov)
            total_inter.append(tin)

            # Calculate strong residues and densities
            strong_resids, _ = self.calculate_strong_resids(trio_pos, utz, ltz, names, resids)
            d2 = self.calculate_densities_for_resids(trio_pos, resids, strong_resids, groups, pbc, bins)
            d2_densities.append(d2)
            sov, sin = self.calculate_overlap_and_inter(d0, d2, dz)
            strong_ov.append(sov)
            strong_inter.append(sin)

            strong_residue_ids.append(strong_resids.tolist()) # convert to list
            with open(f'{self.strong_resid_list_name}.txt', 'w') as f:
                f.write('\n'.join(map(str, strong_residue_ids)))

            

            # Inverted residues densities
            d3 = self.calculate_densities_for_resids(trio_pos, resids, strong_resids, groups, pbc, bins, invert=True)
            d3_densities.append(d3)
            wov, win = self.calculate_overlap_and_inter(d0, d3, dz)
            weak_ov.append(wov)
            weak_inter.append(win)
            strong_num.append(len(strong_resids))
            times.append(ts.time / 1000)
            zs.append(pbc[2] / 10)

        XX = np.linspace(0, np.mean(zs), nbins)
        strong_num = np.array(strong_num)  # Convert to NumPy array for division
        results = {
            'inter': {
                'total': np.transpose([times, total_inter]),
                'strong': np.transpose([times, strong_inter]),
                'weak': np.transpose([times, weak_inter]),
            },
            'ov': {
                'total': np.transpose([XX, np.mean(total_ov, axis=0)]),
                'strong': np.transpose([XX, np.mean(strong_ov, axis=0)]),
                'weak': np.transpose([XX, np.mean(weak_ov, axis=0)]),
            },
            'ratio': {
                'num': np.transpose([times, strong_num]),
                'trio-to-pl': np.transpose([times, strong_num / numP]),
                'trio-to-pl+trio': np.transpose([times, strong_num / (numP + strong_num)]),
            },
            'density': {
                'PL': np.transpose([XX, np.mean(d0_densities, axis=0)]),
                'TRIO': np.transpose([XX, np.mean(d1_densities, axis=0)]),
                'SURF-TRIO': np.transpose([XX, np.mean(d2_densities, axis=0)]),
                'CORE-TRIO': np.transpose([XX, np.mean(d3_densities, axis=0)]),
                'water': np.transpose([XX, np.mean(d_water_series, axis=0)]),
            },
            'strong_residues': strong_residue_ids # List of strong residue IDs
        }

        print("units: Z (nm), interdigitation (nm), time (ns), density (g/m3)")
        return results

    def save_results(self, results, base_dir):
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
        for key1 in results.keys():
            new_dir = os.path.join(base_dir, key1)
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)

            for key2 in results[key1].keys():
                # Define new file path
                file_path = os.path.join(new_dir, f'interdigit_.{key1}.{key2}.dat')
                # Save the file
                np.savetxt(file_path, results[key1][key2])

        print(f"All files have been saved in the directory: {base_dir}")
