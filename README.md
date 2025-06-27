# surface_properties
Under development, code used for Drude paper analysis, soon to be published.

This is used for calculating the surface properties of membrane (bilayer of lipid droplet systems). This includes tail order parameters (Scd), surface triacylglycerol lifetimes, and interdigitation analysis.



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


The example below is for lipid droplet trilayers. However, this can be done with differing bilayer or membrane systems.


![image](https://github.com/user-attachments/assets/e0f088e5-439a-4be5-9205-defefaec8541)


We can then determine the diferent types of interdigition. For example, in a LD monolayer, there is generally triolein interacting with the phospholipids. There would be weak interdigitation (tails of the trioleins interdigitating with the PL), and strong interdigitation (3 or more oxygens of triolein glycerol above z-dim midpoint of PL tails). This is fully modifiable by lipid type, tail chain definition, and strong int. definiton in the functions.

![image](https://github.com/user-attachments/assets/d3f3077e-5357-483c-9cd7-ce5b22383db7)

Strong interdigitation has correlation with the area-per-lipid of monolayer systems:
![image](https://github.com/user-attachments/assets/5162b67b-e262-4baf-8f7f-9cbd95f1b6e4)


# **Lipid Tail Order Parameter (\( S_{cd} \)) Analysis**
The **order parameter** \( S_{cd} \) quantifies the **alignment of lipid tail C-H bonds** relative to the **membrane normal (z-axis)**. It provides insight into:

- **Lipid packing and membrane fluidity**.
- **Bilayer rigidity** under different conditions**.
- **Effects of neutral lipids (e.g., TRIO) on lipid ordering**.

---

## **Order Parameter Definition**
The order parameter is defined as:

$$
S_{cd} = \frac{1}{2} \left( 3 \langle \cos^2 \theta \rangle - 1 \right)
$$

Where:

- \( \theta \) is the **angle between the C-H bond vector and the membrane normal (z-axis)**.
- \( \langle \cdot \rangle \) represents an **ensemble average** over time and molecules.

### **Expected Values**
- \( S_{cd} \approx 1.0 \) → **Highly ordered lipid tails** (rigid packing).
- \( S_{cd} \approx 0.0 \) → **Disordered lipid tails** (fluid-like).
- **Negative values** indicate **tilt away from the normal**.

---

## **Methodology**
This calculation is performed using **MDAnalysis** and follows these steps:

1. **Identify Lipid Tails**: Select **C-H pairs** in lipid acyl chains.
2. **Compute Bond Orientations**: Measure **angles** between **C-H vectors** and the membrane normal.
3. **Average Over Time**: Compute **\( S_{cd} \)** values per **carbon index** across all frames.
4. **Output Results**: Save results in a `.txt` file.

---

## **Example Usage**
### **Running Order Parameter Calculations**
To compute **\( S_{cd} \)** for **POPC tails**, use:

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
T_{\text{surf}} = \sum_{i=0}^{N} \delta_i \cdot \Delta t
$$

where:
- \( T_{\text{surf}} \) is the **total surface residence time**.
- \( \delta_i = 1 \) if the molecule is **surface-bound** at frame \( i \), else **0**.
- \( \Delta t \) is the **frame time step**.
- \( N \) is the **number of trajectory frames**.

---

## **Methodology**
1. **Identify surface molecules** based on \( z \)-position and structural constraints.
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

### **Visualizing Surf-TG Lifetimes**  
```python
import json
import matplotlib.pyplot as plt

with open("lifetime_results/trio_lifetimes.json", "r") as f:
    lifetime_data = json.load(f)

all_lifetimes = [time for lifetimes in lifetime_data.values() for time in lifetimes]

plt.figure(figsize=(6, 4))
plt.hist(all_lifetimes, bins=30, alpha=0.7, color="blue", edgecolor="black")
plt.xlabel("Surface Lifetime (frames)")
plt.ylabel("Frequency")
plt.title("Distribution of Surf-TG Lifetimes")
plt.grid()
plt.show()
```

---

- **Longer lifetimes** indicate **strong surface retention**.
- **Force field-dependent trends** can be observed across simulations.

---
