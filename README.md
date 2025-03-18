# surface_properties
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


![image](https://github.com/user-attachments/assets/e0f088e5-439a-4be5-9205-defefaec8541)
