import logging
import os
import json
logger = logging.getLogger(__name__)

def make_dynamic_residue(strong_residues_list, start_frame=0):
    def selector(ts):
        current_frame = ts.frame - start_frame  # ✅ must come first
        if 0 <= current_frame < len(strong_residues_list):
            residues = strong_residues_list[current_frame]
            logger.info(f"Frame {ts.frame}: Strong residues = {residues}")
            return residues
        else:
            logger.info(f"Frame {ts.frame}: No strong residues available.")
            return []
    return selector


def save_lifetimes_to_json(lifetimes, out_dir, filename='trio_lifetimes.json'):
    if not lifetimes:
        print("No lifetimes to save.")
        return

    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, filename)
    json_data = {str(resid): list(map(int, times)) for resid, times in lifetimes.items()}

    with open(output_path, 'w') as f:
        json.dump(json_data, f)
    print(f"Saved TRIO lifetimes to {output_path}")

