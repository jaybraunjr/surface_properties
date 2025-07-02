import os
import MDAnalysis as mda
from analysis.Surface import InterdigitationAnalysis

ROOT = os.path.dirname(os.path.dirname(__file__))

gro = os.path.join(ROOT, "data", "6.6_2.gro")
xtc = os.path.join(ROOT, "data", "rep1_skip100.xtc")

u = mda.Universe(gro, xtc)

analysis = InterdigitationAnalysis(
    u,
    lipids=["POPC"],
    NL="TRIO",
    water="TIP3",
    start=0,
    stop=10,
)

analysis.run()

print("Total interdigitation:", analysis.results["inter"]["total"][:5])
print("Overlap profile sample:", analysis.results["ov"]["total"][:5])
print("Strong residue count:", analysis.results["ratio"]["num"][:5])
