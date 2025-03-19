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

