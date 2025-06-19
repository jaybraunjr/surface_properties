import logging
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
