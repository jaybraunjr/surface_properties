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
analysis.save_results(out_dir='my_output')

print("Total interdigitation:", analysis.results["inter"]["total"][:5])
print("Overlap profile sample:", analysis.results["ov"]["total"][:5])
print("Strong residue count:", analysis.results["ratio"]["num"][:5])
print("Strong residue count:", analysis.results["density"]["PL"][:5])



# Example 2: DOPE membrane, TRIO neutral lipid
analysis2 = InterdigitationAnalysis(
    u,
    lipids=["DOPE"],
    NL="TRIO",
    water="TIP3",
    tail_atoms=["C316"]
)
analysis2.run()
analysis2.save_results(out_dir='my_output_dope_trio_o11')
print("DOPE/TRIO O11:", analysis2.results["inter"]["total"][:5])


print("Total interdigitation:", analysis2.results["inter"]["total"][:5])
print("Overlap profile sample:", analysis2.results["ov"]["total"][:5])
print("Strong residue count:", analysis2.results["ratio"]["num"][:5])
print("Strong residue count:", analysis2.results["density"]["PL"][:5])
