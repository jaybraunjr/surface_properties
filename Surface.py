from collections import defaultdict
import numpy as np
import os
import json
from tqdm import tqdm

class LifetimeAnalysis:
    def __init__(self, universe, lipids, NL, water, use_ls2=False, use_us2=False, buffer_frames=1, min_oxygens=3, buffer_length=20):
        self.u = universe
        self.lipids = lipids
        self.NL = NL
        self.water = water
        self.use_ls2 = use_ls2
        self.use_us2 = use_us2
        self.buffer_frames = buffer_frames
        self.min_oxygens = min_oxygens  # Minimum number of oxygens required
        self.buffer_length = buffer_length 

    def setup_atom_groups(self):
        halfz = self.u.dimensions[2] / 2
        C2 = ' '.join([f'C2{i}' for i in range(2, 22)])
        C3 = ' '.join([f'C3{i}' for i in range(2, 22)])
        tail_atoms = ' '.join(C2.split() + C3.split())

        lipid_selection = ' or '.join([f'resname {lipid}' for lipid in self.lipids])
        # Define us and ls (standard selections)
        us = f'(same residue as ({lipid_selection}) and name P and prop z>{halfz}) and (name {" ".join(tail_atoms.split())})'
        ls = f'(same residue as ({lipid_selection}) and name P and prop z<{halfz}) and (name {" ".join(tail_atoms.split())})'

        us2 = f'(same residue as ({lipid_selection}) and name C24 and prop z>{halfz})'
        ls2 = f'(same residue as ({lipid_selection}) and name C24 and prop z<{halfz})'

        # Use instance attributes to select upper and lower leaflets
        upper_selection = us2 if self.use_us2 else us
        lower_selection = ls2 if self.use_ls2 else ls

        groups = {
            'memb': self.u.select_atoms(lipid_selection),
            'umemb': self.u.select_atoms(upper_selection),
            'lmemb': self.u.select_atoms(lower_selection),
            'trio': self.u.select_atoms(f'resname {self.NL}'),
            'water': self.u.select_atoms(f'resname {self.water}')
        }
        return groups

    def calculate_strong_resids(self, trio_pos, utz, ltz, names, resids):
        boolArray = ((trio_pos[:, 2] > utz) | (trio_pos[:, 2] < ltz)) & np.char.startswith(names, 'O')
        strong_resids = resids[boolArray]
        r = np.unique(strong_resids, return_counts=True)
        boolArray = (r[1] >= self.min_oxygens) & (r[1] <= 6)
        strong_resids = r[0][boolArray]
        counts = r[1][boolArray]  
        return strong_resids, counts

    def calculate_trio_lifetimes(self, start_frame=0, end_frame=None, step_frame=1):
        if end_frame is None:
            end_frame = self.u.trajectory.n_frames
        groups = self.setup_atom_groups()
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
