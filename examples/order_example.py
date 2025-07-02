import os
import numpy as np
import MDAnalysis as mda
from analysis.order import OrderParameters
import analysis.opc as opc

ROOT = os.path.dirname(os.path.dirname(__file__))

gro = os.path.join(ROOT, "data", "6.6_2.gro")
xtc = os.path.join(ROOT, "data", "rep1_skip100.xtc")

# Load the trajectory
u = mda.Universe(gro, xtc)

# Set up and run the order parameter analysis over the first 10 frames
op = OrderParameters(
    u,
    atomlists=opc.POPC1,
    selection="resname POPC",
    start=0,
    stop=10,
)

op.run()

# Save the carbon index and average S_cd values
out_file = os.path.join(ROOT, "examples", "order_output.txt")
np.savetxt(out_file, op.results["output"], fmt="%d %.6f")
print(f"Saved results to {out_file}")
