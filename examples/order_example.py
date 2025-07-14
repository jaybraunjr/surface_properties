import os
import numpy as np
import MDAnalysis as mda
from analysis.Order import OrderParameters
from analysis.MembraneBase import MembraneAnalysisBase
import analysis.opc as opc
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt


ROOT = os.path.dirname(os.path.dirname(__file__))
gro = os.path.join(ROOT, "data", "6.6_2.gro")
xtc = os.path.join(ROOT, "data", "rep1_skip100.xtc")

u = mda.Universe(gro, xtc)
base = MembraneAnalysisBase(u, ['POPC'], 'TRIO', 'TIP3')  # NL/water are irrelevant here

# Your atomlist for POPC
atomlists = opc.POPC1

# Selections
selection_near = base.select_near_protein("POPC", 10)
selection_donut = base.select_near_protein("POPC", 20, inner_distance=10)
selection_far = base.select_near_protein("POPC", 30, inner_distance=20)
selection_opp_leaflet = f"resname POPC and prop z < {base.halfz}"
selection_bulk = f"resname POPC and prop z > {base.halfz}"

def run_op(selection, label):
    op = OrderParameters(u, atomlists=atomlists, selection=selection, start=0, stop=10)
    op.run()
    print(f"{label} S_cd:", op.results["output"][:, 1])
    return op.results['output']

res_near = run_op(selection_near, "Near Protein")
res_donut = run_op(selection_donut, "Donut (20–10 Å)")
res_far = run_op(selection_far, "Far (30–20 Å)")
res_opp_leaflet = run_op(selection_opp_leaflet, "Opposite Leaflet")
res_bulk = run_op(selection_bulk, "Bulk")

# Plot
plt.figure(figsize=(4, 3))
for res, label in zip(
    [res_near, res_donut, res_far, res_opp_leaflet, res_bulk],
    ["Near Protein", "Donut (20–10 Å)", "Far (30–20 Å)", "Opposite Leaflet", "Bulk"]
):
    plt.plot(res[:, 0], res[:, 1], label=label, lw=2)

plt.xlabel("Carbon Index (Tail Position)")
plt.ylabel("Order Parameter (S)")
plt.title("Order Parameters of POPC under Different Spatial Conditions")
plt.legend()
plt.tight_layout()
plt.savefig("order_parameter_plot.png")
print("Plot saved as order_parameter_plot.png")

