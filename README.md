# ja :curly_loop:
 `UNDER CONSTRUCTION`
## **Introduction**
This utility calculates the classical Jiles-Atherton (JA) magnetic hysteresis model. This model has found wide application in physics and engineering to simulate how ferromagnetic materials become magnetized under the influence of an external field.

The magnetic moment inside a material is affected not only by the external field $(H)$ but also by the magnetization of adjacent regions. This combined field is called the effective field ($H_{e}$):

$$
H_{e} = H + \alpha M
$$

Where $\alpha$ is the domain interaction coefficient.

To describe the "ideal" behavior of a material without hysteresis, the Langevin function ($M_{an}$) is used. It shows what the magnetization would be if there were no defects within the material that impede the movement of the domain walls.

$$
M_{an} = M_{s} \cdot (\coth(\frac{H_{e}}{a}) - \frac{a}{H_{e}})
$$

Real materials contain microscopic defects and impurities. When domain walls (boundaries between magnetic regions) move, they "catch" on these defects. This process is called pinning.
  * To move a wall further, energy must be expended.
  * This expended energy is converted into heat—this is how hysteresis losses occur.
  * In equations, this is described by the coefficient $k$ (pinning factor).

Differentiating the energy balance equation in a magnetic material with respect to $B_{e}$ yields a solution describing the magnetization processes.

$$
M = M_{an} - \delta k \frac{dM}{dB_{e}}
$$

The hysteresis loop can be represented with $B$ as a function of $H$:

$$
B = μ_{0}(H + M)
$$

or with $M$ as a function of $H$, where the dependence of $M$ on $H$ is represented in differential form so that both relationships contain the same information.

The total magnetization ($M$) is divided into two parts:
  * Irreversible ($M_{irr}$): the domain walls overcome obstacles (pinning) and do not return to their original state. This creates a "memory" effect in the material.
  * Reversible ($M_{rev}$): the domain walls simply "bend" under the field pressure, like an elastic membrane, and return to their original state when the field is removed. In the text, this is described by the coefficient $c$.

By combining all factors (energy balance, pinning losses, and elastic deflection of the walls), the authors derive a final differential equation. This allows one to calculate how the magnetization ($dM$) changes with a change in the external field ($dH$).

$$
\frac{dM}{dH} = \frac{1}{1 + c}\frac{1}{\delta k/μ_{0} − \alpha(M_{an} − M)}(M_{an} − M) + \frac{c}{1 + c}\frac{dM_{an}}{dH}
$$

To model any magnetic material (steel, ferrite, nickel), one must substitute into the final equation:

  * $M_{s}$ (Saturation): The maximum possible magnetization. When all the magnetic moments inside are aligned, the material cannot be further magnetized.
  * $a$ (Curve shape): Determines how quickly the material approaches saturation. Depends on the temperature and structure of the material.
  * $k$ (Coercivity / Pinning): Controls the "width" of the loop. The larger $k$, the more defects in the material that the domain walls "catch," and the more energy is lost as heat. 
  * $\alpha$ (Interaction): Indicates how strongly already magnetized regions help (or hinder) the magnetization of their neighbors.
  * $c$ (Reversibility): The "elasticity" coefficient. If $c$ is large, the material behaves like rubber: it easily returns to its shape (magnetic state) when the field is removed.

## **General scheme of the algorithm**
### **Time Preparation and Discretization**
To simulate a periodic signal, a time grid is generated. Unlike simple models, strict indexing is used here:
  * **Signal Period** ($T$): 2.0 seconds.
  * **Simulation Time**: 40 seconds (20 full periods). This is necessary to complete transient processes and allow the model to reach a stable steady-state state.
  * **Integration Step** ($tin$): The recommended value is $0.001$ s to ensure a balance between speed and accuracy.

### **Model Core - Instantaneous Susceptibility**
Instead of the iterative secant method, an analytical calculation of the **hysteresis-free magnetization** ($M_{an}$) and its derivative is used directly within the solver. This is implemented in the `get_dMdH_instant` method, which calculates the slope of the $dM/dH$ curve for any point $(H, M)$.

Key formulas:
  * **Effective field**: $H_e = H + \alpha M$ (takes into account domain interactions).
  * **Langevin function**: Describes the ideal material response. Linearization (Taylor series) is used to protect against singularities at zero ($H \to 0$).
  * **Energy balance**: Irreversible changes ($M_{irr}$) are calculated only under the condition of energy absorption (the condition $(M_{an} - M) \cdot \delta \gt 0$).

### **Runge-Kutta 4th-order integrator (RK4)**
The RK4 method is used to solve the differential equation of magnetization, significantly exceeding the Euler method in accuracy.
At each step, the program calculates four coefficients ($k_1, k_2, k_3, k_4$), "probing" the loop curvature. This allows for an accurate description of the saturation "knee" and minimizes accumulated numerical error.

### **Data Stabilization and Sampling**
The program ignores the first 19 periods ("warm-up") and writes only the last, 20th period to the buffer. This ensures that the resulting loop is steady and symmetrical.

### **Drift Compensation**
Due to microscopic rounding errors in double numbers, the loop may not close perfectly. The algorithm performs the final "stitching":
  * The gap between the first and last point of the period is calculated ($driftH$, $driftM$).
  * The error is distributed linearly across the entire array, making the loop geometrically closed.

### **Final Calculation of Magnetic Flux Density ($B$)**
Magnetic flux density is calculated using the "harmonized" formula from modern research (MDPI):

$$
B = \mu _{0}(H_{e} + M)
$$

Where the magnetization $M$ is pre-centered relative to its average value over the period ($M_{bias}$) to correctly determine the residual flux density ($B_{r}$).

### **Physical Metrics and Analysis**
TODO

### **The flowchart of the `get_dMdH_instant` method**
Describes the logic for calculating the derivative $dM/dH$ to determine the state of a material. This algorithm is an implementation of the J-A physical model in a single time zone.

<img src="https://github.com/aemeltsev/ja/blob/main/img/dmdh_inst.png" width="50%">

### **The flowchart of the RK4 (Runge-Kutta fourth-order) method**

Describes how, within a single time step `tin`, the program performs four "test" derivative calculations
to predict the magnetization trajectory as accurately as possible.

<img src="https://github.com/aemeltsev/ja/blob/main/img/rk4.png" width="50%">

### **The calculate method flowchart**
Describes the high-level logic of the program's operation: from initializing parameters to obtaining a finished, "sewn" hysteresis loop.

<img src="https://github.com/aemeltsev/ja/blob/main/img/calculate.png" width="50%">

## **Recommended reading**
  * Jiye Zhao et al. - State Space Representation of Jiles–Atherton Hysteresis Model and Application (MDPI, 2024).
  * Guangming Xue - Numerical Solving Method for Jiles-Atherton Model and Influence Analysis.
  * Gustav Mörée - Review of Hysteresis Models for Magnetic Materials.
  * Nowicki M. - Modeling the Hysteresis Loop of Ultra-High Permeability Amorphous Alloy.
  * Jastrzebski R. - Comparison of macroscopic descriptions of magnetization curves.
  * Knypinski L. - Application of a PSO algorithm for identification of the parameters of Jiles-Atherton hysteresis model.
  * Kuznetsov V. - Improved Jiles–Atherton Magnetic Core Model and Its SPICE Implementation(MDPI, 2026).
  * Rupnik U. - Harmonization and Validation of Jiles–Atherton Static Hysteresis Models.
  * Szewczyk R.-Validation of the Anhysteretic Magnetization Model for Soft Magnetic Materials with Perpendicular Anisotropy (MDPI, 2014)
  * Gozdur R.-A Study of Temperature-Dependent Hysteresis Curves for a Magnetocaloric Composite Based on La(Fe, Mn, Si)13-H Type Alloys (MDPI, 2020)
