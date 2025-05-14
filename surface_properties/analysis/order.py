
import logging
# from . import opc
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

logger = logging.getLogger(__name__)

class OrderParameters:
    def __init__(self, u, atomlists, selection=None, start_frame=0, end_frame=None, get_strong_residues=None):
        if selection is None and get_strong_residues is None:
            raise ValueError("Either a selection string or a dynamic residue selection callable must be provided.")
        if selection is not None and get_strong_residues is not None:
            raise ValueError("Provide either a selection string or a get_strong_residues callable, not both.")
        if end_frame is not None and end_frame <= start_frame:
            raise ValueError("end_frame must be greater than start_frame")

        self.u = u
        self.atomlists = atomlists
        self.selection = selection
        self.start_frame = start_frame
        self.end_frame = end_frame
        self.get_strong_residues = get_strong_residues
        self.C_numbers, self.Cs, self.Hs_f, self.repeat = self.process_atom_lists()


    def process_atom_lists(self):
        C_numbers = []
        Cs = []
        Hs = []
        repeat = []

        for atoms in self.atomlists:
            # Extract the numeric part of the atom name (e.g., "C1" -> 1)
            C_number = ''.join(filter(str.isdigit, atoms[0]))
            if C_number:
                C_numbers.append(int(C_number))
            else:
                raise ValueError(f"Invalid atom name format: {atoms[0]}")

            Cs.append(atoms[0])
            Hs.append(atoms[1:])
            repeat.append(len(atoms) - 1)

        Hs_f = [item for sublist in Hs for item in sublist]
        assert int(np.sum(repeat)) == len(Hs_f), "Mismatch in repeats"
        return C_numbers, Cs, Hs_f, repeat


    def compute_OP(self):
        if self.end_frame is not None and self.end_frame <= self.start_frame:
            raise ValueError("end_frame must be greater than start_frame")
        output = []
        total_frames = len(self.u.trajectory) if self.end_frame is None else (self.end_frame - self.start_frame)
        for idx, ts in enumerate(self.u.trajectory[self.start_frame:self.end_frame]):
            logger.info(f"Processing frame {idx+1} of {total_frames}")
            if self.get_strong_residues is not None:
                strong_residues = self.get_strong_residues(ts)
                if len(strong_residues) == 0:
                    output.append(np.zeros(len(self.C_numbers)))
                    continue
                strong_res_selection = "resid " + " or resid ".join(map(str, strong_residues))
                all_molecules = self.u.select_atoms(strong_res_selection)
            else:
                all_molecules = self.u.select_atoms(self.selection)
            if len(all_molecules) == 0:
                logger.info(f"No atoms found for frame {idx+1}. Skipping.")
                output.append(np.zeros(len(self.C_numbers)))
                continue
            valid_indices_group1 = []
            valid_indices_group2 = []
            for molecule in all_molecules.residues:
                group1_atoms = molecule.atoms.select_atoms("name " + " ".join(self.Cs))
                group2_atoms = molecule.atoms.select_atoms("name " + " ".join(self.Hs_f))
                if len(group1_atoms) > 0 and len(group2_atoms) > 0:
                    valid_indices_group1.extend(group1_atoms.indices)
                    valid_indices_group2.extend(group2_atoms.indices)
            if len(valid_indices_group1) == 0 or len(valid_indices_group2) == 0:
                output.append(np.zeros(len(self.C_numbers)))
                continue
            group1 = self.u.atoms[valid_indices_group1]
            group2 = self.u.atoms[valid_indices_group2]
            natoms = len(self.Cs)
            nmols = int(len(group1.positions) / natoms)
            repeats = self.repeat * nmols
            p1 = np.repeat(group1.positions, repeats, axis=0)
            p2 = group2.positions
            dp = p2 - p1
            norm = np.sqrt(np.sum(np.power(dp, 2), axis=-1))
            norm[norm == 0] = 1e-8
            cos_theta = dp[..., 2] / norm
            S = -0.5 * (3 * np.square(cos_theta) - 1)
            new_S = self._average_over_hydrogens(S, repeats)
            new_S.shape = (nmols, natoms)
            results = np.average(new_S, axis=0)
            output.append(results)
        if len(output) == 0:
            logger.info("No results.")
            return np.array([])
        avg = np.average(output, axis=0)
        return np.transpose([self.C_numbers, avg])

    def _average_over_hydrogens(self, x, reps):
        i = 0
        out = []
        for rep in reps:
            tmp = x[i:i+rep]
            i += rep
            out.append(np.average(tmp))
        return np.array(out)

def run_op(u, opc, lipid_selection, selection=None, start_frame=0, end_frame=None, output_text='pl', get_strong_residues=None):
    if not hasattr(opc, lipid_selection):
        raise ValueError(f"No lipid atoms configuration found for {lipid_selection}")
    lipid_atoms = getattr(opc, lipid_selection)
    OP_calc = OrderParameters(u, lipid_atoms, selection=selection, start_frame=start_frame, end_frame=end_frame, get_strong_residues=get_strong_residues)
    OP_results = OP_calc.compute_OP()
    if OP_results.size > 0:
        np.savetxt(output_text, OP_results)
    else:
        logger.info("No order parameter results to save.")
    return OP_results


def run_vec_leaflet(
    universe,
    frame_start,
    frame_end,
    frame_step,
    residue_sel="resname TRIO",
    tail_names=['C118', 'C218', 'C318'],
    headgroup_sel="name O*",
    pl_selection="resname POPC DOPE and name C210",
    leaflet='bottom',  # 'bottom' = z > midpoint, 'top' = z < midpoint
    expected_headgroup_count=6
):
    angles = []
    time_series = []
    avg_order_params = {tail: [] for tail in tail_names}
    std_order_params = {tail: [] for tail in tail_names}

    total_frames = len(universe.trajectory[frame_start:frame_end:frame_step])
    for ts in tqdm(universe.trajectory[frame_start:frame_end:frame_step], total=total_frames,
                   desc=f"Analyzing {universe.trajectory.filename.split('/')[-1]}"):

        pl_atoms = universe.select_atoms(pl_selection)
        if len(pl_atoms) == 0:
            continue

        z_coords = pl_atoms.positions[:, 2]
        midpoint_z = np.mean(z_coords)
        if leaflet == 'bottom':
            leaflet_atoms = pl_atoms.select_atoms(f'prop z > {midpoint_z}')
        elif leaflet == 'top':
            leaflet_atoms = pl_atoms.select_atoms(f'prop z < {midpoint_z}')
        else:
            raise ValueError("leaflet must be 'bottom' or 'top'")
        if len(leaflet_atoms) == 0:
            continue
        avg_leaflet_z = np.mean(leaflet_atoms.positions[:, 2])
        residues = universe.select_atoms(residue_sel).residues
        frame_cosine_values = {tail: [] for tail in tail_names}

        for res in residues:
            headgroup = res.atoms.select_atoms(headgroup_sel)
            if len(headgroup) != expected_headgroup_count:
                continue
            com_head = headgroup.center_of_mass()

            if (leaflet == 'bottom' and com_head[2] <= avg_leaflet_z) or \
               (leaflet == 'top' and com_head[2] >= avg_leaflet_z):
                continue

            for tail in tail_names:
                tail_atom = res.atoms.select_atoms(f'name {tail}')
                if len(tail_atom) != 1:
                    continue
                tail_pos = tail_atom.positions[0]
                vec = com_head - tail_pos
                norm = np.linalg.norm(vec)
                if norm == 0:
                    continue
                cos_theta = np.dot(vec, [0, 0, 1]) / norm

                angles.append({
                    'frame': ts.frame,
                    'resid': res.resid,
                    'tail_name': tail,
                    'cosine_alignment': cos_theta
                })
                frame_cosine_values[tail].append(cos_theta)

        for tail in tail_names:
            values = frame_cosine_values[tail]
            if values:
                avg_order_params[tail].append(np.mean(values))
                std_order_params[tail].append(np.std(values))
            else:
                avg_order_params[tail].append(np.nan)
                std_order_params[tail].append(np.nan)

        time_series.append(ts.frame)

    return angles, time_series, avg_order_params, std_order_params
