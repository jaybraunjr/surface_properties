# **surface\_properties**

This package calculates **surface properties of membranes**, with a focus on **lipid droplet systems**, **bilayers**, and **monolayers**. It includes modular, extensible classes for analyzing **interdigitation**, **tail ordering**, **cosine alignment**, and **surface molecule lifetimes**. It is designed for use with **MDAnalysis** and supports customization of lipid types, molecular definitions, and spatial criteria.

---
## **Core Analyses**

### **1. Interdigitation Analysis**

Quantifies the overlap between **neutral lipids** (e.g. TRIO) and **phospholipids** based on **z-direction density profiles**.

* Computes **total**, **strong**, and **weak interdigitation** based on glycerol oxygen positioning.
* Defines overlap parameter $\rho_{ov}(z) \in [0, 1]$, and integrates over $z$ to get interdigitation $\lambda_{ov}$.
* Supports customizable leaflet definitions, tail atoms, and binning resolution.
* Modular design using `InterdigitationAnalysis` inheriting from `AnalysisBase`.

### **2. Tail Order Parameter $S_{cd}$**

Measures **lipid tail ordering** with respect to the membrane normal (z-axis).

* Calculates $S_{cd} = \frac{1}{2}(3\langle \cos^2\theta \rangle - 1)$ for each carbon position.
* Allows selections near/away from protein, or by residue.
* Easily extended for any hydrocarbon tail type.

### **3. Cosine Alignment Analysis**

Tracks the **time-dependent orientation** of tail vectors (e.g. C–H or terminal–mid tail vectors) relative to z-axis.

* Outputs per-frame cosine values and averages.
* Suitable for autocorrelation and alignment decay analysis.

### **4. Surface Molecule Lifetime**

Calculates **surface residence times** of neutral lipids (e.g. surface-active TGs).

* Tracks per-residue surface binding using z-position and buffer thresholds.
* Returns lifetimes per molecule and allows histogramming or export to `.json`.

---

## **Modularity & Workflow**

All analyses inherit from a unified `AnalysisBase` class:

* Shared `run()` interface for consistency.
* Flexible per-frame `_analyze_frame()` logic.

You can run all analyses in just a few lines of code:

```python
from analysis.order import OrderParameters

op = OrderParameters(
    universe=u,
    atomlists=opc.POPC1,                        # list of CH bonds
    selection="resname POPC and prop z > 6.0",  # e.g., top leaflet
    start=0,
    stop=100
)

op.run()
results = op.results['output']  # Columns: [carbon_index, Scd]


from analysis.Vector import VectorOrientation

vec = VectorOrientation(
    u,
    residue_sel="resname TRIO",
    tail_names=["C118", "C218", "C318"],
    pl_selection="resname POPC and name C210",
    leaflet="bottom",
    start=0,
    stop=100
)

vec.run()
angles, time_series, avg_order, std_order = vec.unpack()



from analysis.Surface import InterdigitationAnalysis

inter = InterdigitationAnalysis(
    u,
    lipids=['POPC'],
    NL='TRIO',
    water='TIP3',
    start=0,
    stop=100
)

results = inter.run().results


from analysis.lifetime import LifetimeAnalysis

lt = LifetimeAnalysis(
    universe=u,
    lipids=["POPC"],
    NL="TRIO",
    water="TIP3",
    buffer_frames=3,
    min_oxygens=3,
    buffer_length=20
)

lifetimes = lt.calculate_trio_lifetimes(start_frame=0, end_frame=100)
lt.analyze_and_save(base_dir="results/")

```

---

## **Customization**

* Supports **custom tail atom definitions**.
* Compatible with **trilayers**, **bilayers**, and **monolayer geometries**.

