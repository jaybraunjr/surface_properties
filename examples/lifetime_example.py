import os
import MDAnalysis as mda
from analysis.Lifetime import LifetimeAnalysis
from analysis.utils import save_lifetimes_to_json

ROOT = os.path.dirname(os.path.dirname(__file__))

gro = os.path.join(ROOT, "data", "6.6_2.gro")
xtc = os.path.join(ROOT, "data", "rep1_skip100.xtc")

u = mda.Universe(gro, xtc)

analysis = LifetimeAnalysis(
    u,
    lipids=["POPC"],
    NL="TRIO",
    water="TIP3",
    start=0,
    stop=10,
)

analysis.run()

out_dir = os.path.join(ROOT, "examples", "lifetime_results")
save_lifetimes_to_json(analysis.results, out_dir)
