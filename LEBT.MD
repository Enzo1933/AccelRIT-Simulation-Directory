# RITACCEL Theory Manual
## Low Energy Beam Transport (LEBT) & Electrostatic Optics

### Introduction
The Low Energy Beam Transport (LEBT) module of RITACCEL is designed to handle highly divergent, non-relativistic ion beams extracted from an ion source (e.g., $20\text{ keV}$ $H^-$ ions). The primary focusing element in this regime is the **Einzel Lens**, a symmetric three-electrode electrostatic device. By raising the central electrode to a high potential while grounding the outer electrodes, the lens focuses the beam radially without altering its final kinetic energy.

### The Analytical Field Model
To avoid the computational overhead of importing Finite Element Analysis (FEA) meshes for routine optimization, RITACCEL uses a continuous analytical approximation for the electric potential of a symmetric three-cylinder lens.

For a lens with cylinder radius $R$ and a middle electrode of length $L$ sitting at voltage $U$, we define a geometric constant $k = 1.318 / R$. The front and rear electrode gaps are defined longitudinally as:
$$u_1 = k\left(z + \frac{L}{2}\right)$$
$$u_2 = k\left(z - \frac{L}{2}\right)$$ 

The on-axis potential $\Phi(z)$ and the longitudinal electric field $E_z$ (which is $-\Phi'(z)$) are evaluated as:
$$\Phi(z) = \frac{U}{2} \left[ \tanh(u_1) - \tanh(u_2) \right]$$
$$\Phi'(z) = \frac{U \cdot k}{2} \left[ \text{sech}^2(u_1) - \text{sech}^2(u_2) \right]$$

Because the non-relativistic beam decelerates (or accelerates) inside the lens, the total kinetic energy equivalent voltage $V(z)$ of the beam at any point is strictly tracked as:
$$V(z) = V_0 + \Phi(z)$$
*(Where $V_0$ is the initial extraction voltage of the ion source).*

### The Paraxial Ray Equation & The Picht Substitution
The trajectory of a particle in an electrostatic field with rotational symmetry is governed by the Paraxial Ray Equation:
$$r'' + \frac{V'(z)}{2V(z)}r' + \frac{V''(z)}{4V(z)}r = 0$$

Numerically integrating this equation using standard methods (like Runge-Kutta) is highly sensitive to step size because of the first-derivative term ($r'$), which represents the beam's changing axial velocity. 

To achieve absolute numerical stability, RITACCEL applies the **Picht Substitution** (Quadrature), transforming the physical radius $r(z)$ into a "reduced" coordinate $R(z)$:
$$R(z) = r(z) \cdot [V(z)]^{1/4}$$

Substituting $R(z)$ back into the paraxial equation perfectly eliminates the $r'$ velocity term, yielding the **Reduced Paraxial Ray Equation**:
$$R'' = -\frac{3}{16} \left(\frac{V'(z)}{V(z)}\right)^2 R$$

**Physics Implication:** Because the coefficient on the right side of the equation is strictly negative, the resulting force is *always* a restoring force toward the center axis. This provides the mathematical guarantee that an Einzel lens is always focusing, regardless of the applied voltage polarity.

### The Numerical Engine: Piecewise Matrix Slicing
Notice that the Reduced Paraxial Ray Equation takes the exact form of a **Simple Harmonic Oscillator**: $x'' = -\omega^2 x$. The "focusing frequency" at any localized point $z$ is:
$$\omega(z) = \frac{\sqrt{3}}{4} \left| \frac{V'(z)}{V(z)} \right|$$

Because $\omega$ changes continuously, RITACCEL does not attempt a global analytical solution. Instead, it utilizes **Piecewise Constant Slicing**. The engine chops the $Z$-axis into microscopic slices ($dz$). Inside each $dz$ slice, $\omega$ is assumed to be momentarily constant, allowing the engine to propagate the beam state using exact, unconditionally stable Transfer Matrices:

$$\begin{pmatrix} R_{next} \\ R'_{next} \end{pmatrix} = \begin{pmatrix} \cos(\omega \, dz) & \frac{1}{\omega} \sin(\omega \, dz) \\ -\omega \sin(\omega \, dz) & \cos(\omega \, dz) \end{pmatrix} \begin{pmatrix} R_{current} \\ R'_{current} \end{pmatrix}$$

(Note: If the local electric field V'(z) approaches zero, $\omega \to 0$, and the engine dynamically collapses the transport matrix into a standard vacuum drift matrix).

At the end of the simulation pipeline, the final reduced coordinates $R_f$ and $R'_f$ are converted back into the physical radius $r_f$ and divergence angle $r'_f$ for accurate phase-space handoff to downstream acceleration columns.