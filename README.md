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

Got it! Here’s your **README documentation** in the **correct LaTeX Markdown format** for direct pasting:

---

# **Lipid Tail Order Parameter (\( S_{cd} \)) Analysis**
The **order parameter \( S_{cd} \)** quantifies the **alignment of lipid tail C-H bonds** relative to the **membrane normal (z-axis)**. It provides insight into:
- **Lipid packing and membrane fluidity**.
- **Bilayer rigidity** under different conditions.
- **Effects of neutral lipids (e.g., TRIO) on lipid ordering**.

## **Order Parameter Definition**
The order parameter is defined as:

$$
S_{cd} = \frac{1}{2} \left( 3 \langle \cos^2 \theta \rangle - 1 \right)
$$

Where:
- \( \theta \) is the **angle between the C-H bond vector and the membrane normal (z-axis)**.
- \( \langle \cdot \rangle \) represents an **ensemble average** over time and molecules.

### **Expected Values:**
- \( S_{cd} \approx 1.0 \) → **Highly ordered lipid tails** (rigid packing).
- \( S_{cd} \approx 0.0 \) → **Disordered lipid tails** (fluid-like).
- **Negative values** indicate tilt **away from the normal**.

---

## **Methodology**
This calculation is performed using **MDAnalysis** and follows these steps:

1. **Identify Lipid Tails:** Select **C-H pairs** in lipid acyl chains.
2. **Compute Bond Orientations:** Measure **angles** between **C-H vectors** and the membrane normal.
3. **Average Over Time:** Compute **\( S_{cd} \)** values per **carbon index** across all frames.
4. **Output Results:** Save results in a `.txt` file.

---

## **Example Usage**
### **Running Order Parameter Calculations**
To compute **\( S_{cd} \)** for **POPC tails**, use:

```python
import MDAnalysis as mda
from order import run_op

# Load trajectory
universe = mda.Universe("membrane.gro", "trajectory.xtc")

# Run order parameter calculation for POPC tails
order_results = run_op(
    u=universe,
    opc=opc,  
    lipid_selection="POPC",  # Define lipid of interest
    selection="resname POPC",
    start_frame=0,
    end_frame=500,
    output_text="order_parameters.txt"
)
```
This will generate a text file **`order_parameters.txt`** containing **\( S_{cd} \) values per carbon index**.

---

## **Interpreting the Results**
The output contains **two columns**:
- **Column 1:** Carbon index \( C_i \).
- **Column 2:** Order parameter \( S_{cd} \).

Example:
```
1   0.23
2   0.35
3   0.52
...
```

### **Visualizing the Order Parameter Profile**
You can **plot** \( S_{cd} \) profiles using:

```python
import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt("order_parameters.txt")
carbon_numbers, scd_values = data[:, 0], data[:, 1]

# Plot
plt.figure(figsize=(6, 4))
plt.plot(carbon_numbers, scd_values, marker='o', linestyle='-', label="POPC")
plt.xlabel("Carbon Index")
plt.ylabel(r"$S_{cd}$")
plt.title("Lipid Tail Order Parameter")
plt.legend()
plt.grid()
plt.show()
```

### **What to Look For:**
- **High \( S_{cd} \) (~0.8–1.0) at tail ends** → Tighter packing.
- **Lower \( S_{cd} \) at mid-chain** → More disorder.
- **Drop in \( S_{cd} \) near the headgroup** → Increased mobility.

---

## **Why Is This Important?**
- **TRIO insertion into lipid monolayers** **reduces \( S_{cd} \)** → increased fluidity.
- **Bilayers with high \( S_{cd} \) are more rigid**, affecting lipid mobility.
- **Different lipids (e.g., POPC vs. DOPE)** show distinct **\( S_{cd} \)** profiles, impacting biological function.

---

## **Future Extensions**
- Compute **\( S_{cd} \)** for different lipid species (e.g., DOPE, DPPC).
- Compare **TRIO-exposed vs. TRIO-free regions**.
- Analyze **bilayers vs. monolayers**.

---

This is now in the **exact format** you can **paste directly into your README** with proper **Markdown LaTeX formatting**! 🚀  

Would you like me to add **any specific examples** or **references to lipid order in literature**?

