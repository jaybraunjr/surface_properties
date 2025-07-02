import os
import MDAnalysis as mda
from analysis.Vector import VectorOrientation

ROOT = os.path.dirname(os.path.dirname(__file__))

gro = os.path.join(ROOT, "data", "6.6_2.gro")
xtc = os.path.join(ROOT, "data", "rep1_skip100.xtc")

u = mda.Universe(gro, xtc)

vo = VectorOrientation(u, start=0, stop=10)
vo.run()

angles, times, avg_order, std_order = vo.unpack()
print("First angle entry:", angles[0])
print("Time series length:", len(times))
print("Avg order keys:", list(avg_order.keys()))
