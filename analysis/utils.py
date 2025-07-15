"""Utility helpers for saving and selecting residues.

The helper functions in this module assist with common tasks such as
dynamically selecting residue ids per frame and serialising results to
JSON files.
"""

import logging
import os
import json
logger = logging.getLogger(__name__)

def make_dynamic_residue(strong_residues_list, start_frame=0):
    """Return a selector that yields strong residues for each frame.

    Parameters
    ----------
    strong_residues_list : list[list[int]]
        Lists of residue ids deemed strong for each frame.
    start_frame : int, optional
        Frame number corresponding to ``strong_residues_list[0]``.

    Returns
    -------
    Callable
        Function that accepts a :class:`Timestep` and returns a list of residue
        ids.
    """

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
    """Persist lifetimes dictionary to a JSON file.

    Parameters
    ----------
    lifetimes : dict[int, list[int]]
        Mapping of residue id to a list of lifetime lengths.
    out_dir : str
        Directory where the JSON file will be written.
    filename : str, optional
        Name of the output file.
    """

    if not lifetimes:
        print("No lifetimes to save.")
        return

    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, filename)
    json_data = {str(resid): list(map(int, times)) for resid, times in lifetimes.items()}

    with open(output_path, 'w') as f:
        json.dump(json_data, f)
    print(f"Saved TRIO lifetimes to {output_path}")

