# surface_properties


This is used for calculating the surface properties of membrane (bilayer/monolayer or lipid droplet systems). This includes interdigitation analysis, tail order parameters (Scd), surface molecule lifetimes.

# Interdigitation Analysis

The overlap parameter, $\rho_{ov}(z)$, ranges from **0 to 1**, where:
- **0** indicates no overlap.
- **1** represents identical density distributions at a given $z$-position.

It is defined as:

$$
\rho_{ov}(z) = 4 \frac{\rho_{TG}(z) \rho_{PL}(z)}{(\rho_{TG}(z) + \rho_{PL}(z))^2}
$$

The **total interdigitation** $\lambda_{ov}$ is obtained by integrating the overlap parameter along the $z$-axis:

$$
\lambda_{ov} = \int_0^L \rho_{ov}(z) \,dz
$$

## **Methodology**

The interdigitation analysis quantifies how lipid tails and surface molecules (e.g., triolein) interpenetrate along the membrane normal (z-axis). It follows these steps:

1. **Define Groups**  
   Select **phospholipid tails** and **neutral lipid atoms** (e.g., TRIO). Optionally, define *strong interdigitation* criteria, such as a minimum number of oxygen atoms above the PL tail midplane.
2. **Compute Density Profiles**  
   Calculate normalized **z-density distributions** for both **PL** and **TG** groups.
3. **Calculate Overlap Parameter**  
   At each z-position, compute the **overlap** ρ<sub>ov</sub>(z)
4. **Integrate to Get Total Interdigitation**  
Integrate ρ<sub>ov</sub>(z) across the z-dimension to yield:


The example below is for lipid droplet (LD) trilayers. However, this can be done with differing bilayer or membrane systems.

![image](https://github.com/user-attachments/assets/e0f088e5-439a-4be5-9205-defefaec8541)

We can then determine the diferent types of interdigition. For example, in a LD monolayer, there is generally triolein interacting with the phospholipids. There would be weak interdigitation (tails of the trioleins interdigitating with the PL), and strong interdigitation (3 or more oxygens of triolein glycerol above z-dim midpoint of PL tails). This is fully modifiable by lipid type, tail chain definition, and strong int. definiton in the functions.


<table>
<tr>
<td>
  <img src="https://github.com/user-attachments/assets/5162b67b-e262-4baf-8f7f-9cbd95f1b6e4" width="180"/>
</td>
<td>
  <b>Strong interdigitation</b> has implications to biological systems. For example, it correlates with the area-per-lipid of membranes, leading to larger membrane defects.
</td>
</tr>
</table>

## **Example Usage**
To compute **lipid interdigitation**, use:

```python
import MDAnalysis as mda
from analysis.Surface import InterdigitationAnalysis

# Load trajectory
u = mda.Universe("membrane.gro", "trajectory.xtc")

# Set up and run the analysis
analysis = InterdigitationAnalysis(
    u,
    lipids=["POPC"],        # Membrane lipid(s)
    NL="TRIO",              # Neutral lipid residue name
    water="TIP3",           # Water residue name
    start=0,
    stop=500,               # Analyze frames 0–500
    tail_atoms=None,        # Use default membrane tail atoms
    min_oxygens=3           # Strong interdigitation: 3+ oxygens
)
analysis.run()
analysis.save_results(out_dir="interdig_results")

# Access results in Python
print("Total interdigitation:", analysis.results["inter"]["total"][:5])
print("Overlap profile sample:", analysis.results["ov"]["total"][:5])
print("Strong/weak counts:", analysis.results["ratio"][f"{analysis.base.NL.lower()}-to-pl"][:5])
```


# Lipid Tail Order Parameter ($S_{cd}$) Analysis

The **order parameter** $S_{cd}$ quantifies the **alignment of lipid tail C-H bonds** relative to the **membrane normal (z-axis)**.

---

## **Order Parameter Definition**
The order parameter is defined as:

$$
S_{cd} = \frac{1}{2} \left( 3 \langle \cos^2 \theta \rangle - 1 \right)
$$

Where:

- θ is the **angle between the C–H bond vector and the membrane normal (z-axis)**.
- ⟨·⟩ represents an **ensemble average** over time and molecules.

### **Expected Values**
- S<sub>cd</sub> ≈ 1.0 → **Highly ordered lipid tails** (rigid packing).
- S<sub>cd</sub> ≈ 0.0 → **Disordered lipid tails** (fluid-like).
- **Negative values** indicate **tilt away from the normal**.

---
## **Methodology**
This calculation is performed using **MDAnalysis** and follows these steps:

1. **Identify Lipid Tails**: Select **C–H pairs** in lipid acyl chains.
2. **Compute Bond Orientations**: Measure **angles** between **C–H vectors** and the membrane normal.
3. **Average Over Time**: Compute **S<sub>cd</sub>** values per **carbon index** across all frames.
4. **Output Results**: Save results in a `.txt` file.

---

## **Example Usage**
### **Running Order Parameter Calculations**
To compute **S<sub>cd</sub>** for **POPC tails**, use:

```python
import MDAnalysis as mda
import numpy as np
from analysis.order import OrderParameters
from analysis import opc

# Load trajectory
u = mda.Universe("membrane.gro", "trajectory.xtc")

# Set up and run the analysis
op = OrderParameters(
    u,
    atomlists=opc.POPC1,
    selection="resname POPC",
    start=0,
    stop=500
)
op.run()
np.savetxt("order_parameters.txt", op.results["output"], fmt="%d %.6f")
```
This generates **`order_parameters.txt`** with **\( S_{cd} \)** values per carbon index.

---

Here’s the **concise, LaTeX-formatted** README section for **Surf-TG Lifetime Analysis**, ready for direct pasting:

---

# **Surface-Active Molecule Lifetime Analysis**  

The **lifetime analysis** tracks **surface-active molecules** at the **lipid interface**. While focused on **triolein (Surf-TG)**, this method applies to **any molecule** interacting with **monolayers, bilayers, or membranes**.  

### **Defining Surface Persistence**  
A molecule is **surface-active** if it remains at the interface for **continuous frames**. The total surface residence time is:

$$
T_{\text{surf}} = \sum_{i=0}^{N-1} \delta_i \Delta t
$$

where:
- T<sub>surf</sub> is the **total surface residence time**.
- δ<sub>i</sub> = 1 if the molecule is **surface-bound** at frame <i>i</i>, else **0**.
- Δt is the **frame time step**.
- N is the **number of trajectory frames**.

---

## **Methodology**
1. **Identify surface molecules** based on z-position and structural constraints.
2. **Track residence states** across trajectory frames.
3. **Apply a buffer period** to remove transient fluctuations.
4. **Store residence times** for each molecule.

---

## **Example Usage**
To compute **Surf-TG lifetimes**, use:

```python
import MDAnalysis as mda
from analysis.Lifetime import LifetimeAnalysis
from analysis.utils import save_lifetimes_to_json

# Load trajectory
universe = mda.Universe("membrane.gro", "trajectory.xtc")

# Initialize lifetime analysis
lifetime_analysis = LifetimeAnalysis(
    universe,
    lipids=["POPC"],
    NL="TRIO",
    water="TIP3",
    buffer_length=20,
    min_oxygens=3
)

# Run the analysis over frames 0–500
lifetime_analysis.run()

# Save results to JSON
save_lifetimes_to_json(lifetime_analysis.results, "lifetime_results")
```

This creates **`lifetime_results/trio_lifetimes.json`**, storing lifetimes per molecule.

---

## **Interpreting the Results**  
The output contains **surface residence times per molecule**:

```json
{
    "101": [10, 12, 18],
    "205": [5, 7, 25],
    "312": [30, 40, 42]
}
```

where keys represent **molecule IDs** and values are **lifetimes**.

---

- **Longer lifetimes** indicate **strong surface retention**.
- **Force field-dependent trends** can be observed across simulations.

---

## Examples

The repository includes example scripts demonstrating the analysis classes.
Run them from the repository root with Python:

```bash
python examples/order_example.py
python examples/lifetime_example.py
python examples/interdigitation_example.py
python examples/vector_orientation_example.py
```

Each script loads the data in `data/` and prints or saves results in the
`examples/` folder.
