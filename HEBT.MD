# RITACCEL Theory Manual
## High Energy Beam Transport (HEBT) & Magnetic Optics

### Introduction
The High Energy Beam Transport (HEBT) module simulates the beamline after the primary acceleration column. In this regime (e.g., $1\text{ MeV}+$ protons), the beam is stiff enough that electrostatic focusing becomes impractical due to the massive voltages required. Instead, RITACCEL utilizes **Magnetic Quadrupole Triplets** to achieve sub-micron point focusing on the final target.

Unlike the LEBT module, which tracks changing beam velocities, the HEBT module assumes the beam has reached its final plateau energy. The focus shifts to handling highly non-linear, 3D real-world magnetic constraints.

### The Physics: The Frohlich-Kennelly Saturation Model
Standard matrix-based beam optics codes (like basic TRACE3D implementations) assume a linear relationship between the applied coil current (Ampere-turns, $NI$) and the magnetic gradient ($G$). This allows purely mathematical optimizers to request impossible gradients (e.g., $150\text{ T/m}$), which would require the iron pole tips to reach non-physical saturation levels ($B > 2.5\text{ T}$).

RITACCEL prevents "mathematical ghosts" by enforcing the **Frohlich-Kennelly (F-K) Saturation Model**. The magnetic flux density $B$ in the iron yoke is strictly bounded:
$$B(H) = \frac{\mu_0 H}{1 + c H}$$
Where $c$ is a saturation parameter specific to the yoke material (e.g., 1010 Steel). 

By integrating the F-K curve into the core geometry solver, RITACCEL dynamically calculates the true, saturated gradient $G_{true}$ for any given MMF ($NI$). If the numerical optimizer requests too much current, the iron saturates, the gradient flatlines, and the solver mathematically recognizes a "dead end," guaranteeing that all optimized outputs are physically buildable.

### The Physics: Enge Functions & Fringe Fields
A standard "hard-edge" magnet model perfectly contains the magnetic field within the physical length of the iron yoke. In reality, magnetic flux leaks out into the beam pipe, meaning the beam begins focusing *before* it enters the magnet and continues focusing *after* it leaves.

To achieve precision sub-micron tracking, RITACCEL models the smooth spatial decay of the gradient $G(z)$ using **Enge Functions**:
$$F(z) = \frac{1}{1 + \exp\left(a_1 + a_2\left(\frac{z}{D}\right) + a_3\left(\frac{z}{D}\right)^2 + ...\right)}$$
*(Where $D$ is the full aperture of the magnet, and $a_i$ are empirical fall-off coefficients).*

Combined with the **Roters Fringing Factor**, the engine accurately calculates the effective magnetic length ($L_{eff}$), significantly shifting the predicted focal point compared to idealized models.

### The Mathematical Engine: 4th-Order Runge-Kutta (RK4)
Because the magnetic gradient $G(z)$ is now a highly non-linear, continuous curve (dictated by Enge fringe fields and F-K saturation), Transfer Matrices are no longer sufficient. 

The transverse motion of the particles in the X and Y planes is governed by Hill's Differential Equations:
$$x'' + k_x(z)x = 0$$
&
$$y'' - k_y(z)y = 0$$
Where the focusing strength $k(z) = \frac{q \cdot G(z)}{p}$ depends on the particle's charge $q$ and momentum $p$.

RITACCEL integrates these coupled differential equations using a custom **4th-Order Runge-Kutta (RK4)** solver. The RK4 integrator evaluates the derivatives at four distinct points within every microscopic step $dz$, ensuring absolute spatial fidelity as the beam rides down the non-linear fringe field curves.

### The Optimizer: Trust-Region Newton-Raphson
The core feature of the HEBT module is its ability to autonomously find the exact Ampere-turns ($MMF_1, MMF_2$) required for a Quadrupole Triplet to focus the beam perfectly on target. 

To solve this inverse problem, RITACCEL deploys a **Newton-Raphson Root-Finding Algorithm** aimed at minimizing a residual vector (the final $X$ and $Y$ beam sizes at the target distance). 

Because the F-K saturation curve creates "walls" in the parameter space, a standard Newton solver would overshoot, calculate a Singular Jacobian, and crash. RITACCEL mitigates this by applying a **Trust-Region Speed Limit**. The optimizer is mathematically prohibited from changing the MMF by more than 20% per iteration. This forces the solver to gently "walk down the valley" of the gradient field, resulting in a rock-solid, crash-proof convergence to $10^{-8}$ tolerances.