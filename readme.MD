# RITACCEL ⚛️
**Accelerator Computing & Control Engine Library**

RITACCEL is a high-performance, physics-aware particle tracking and optimization engine built in Rust. It is designed to simulate continuous beamlines for RIT's AccelRIT Tandem Accelerator, with a team led by Enzo Adorno.

RITACCEL enforces "hard physical constraints," ensuring that numerical optimizers cannot recommend physically impossible magnet gradients or perfectly sharp electrostatic fields.

## Core Capabilities

The entire engine is split into highly specialized modules, united by a continuous `Beam` state tracker. **As of now, the engine handles 2/4 intended regions of the particle accelerator.**

### 1. Low Energy Beam Transport (LEBT) 
Handles highly divergent, low-energy ion beams (e.g., $20\text{ keV}$ $H^-$) using electrostatic Einzel lenses.
* **Physics:** Uses continuous analytical voltage models ($\tanh$/$\text{sech}$) rather than hard-edge approximations.
* **Math Engine:** Applies the Picht Substitution to the Paraxial Ray Equation, allowing the engine to track changing beam velocities unconditionally.
* **Integrator:** Lightning-fast **Piecewise Transfer Matrix Slicer** simulating localized Simple Harmonic Oscillator optics.

### 2. High Energy Beam Transport (HEBT) 
Optimizes focusing pipelines for high-energy beams using Magnetic Quadrupole Triplets.
* **Physics:** Enforces **Frohlich-Kennelly Magnetic Saturation**, guaranteeing the optimizer cannot request infinite iron permeability. Includes **Enge Functions** and the **Roters Fringing Factor** for realistic 3D flux leakage.
* **Math Engine:** **4th-Order Runge-Kutta (RK4)** integration for sub-micron spatial tracking.
* **Optimizer:** Custom Trust-Region Newton-Raphson solver utilizing analytical Jacobians to find exact Ampere-turn (MMF) requirements for a perfect focal spot.

## 🧠 The Architecture Pipeline
RITACCEL is designed to pass a single particle state through totally distinct physics environments without losing precision. 

## 🛣️ Development Roadmap
- [x] RK4 Magnetic Tracking
- [x] Trust-Region MMF Optimization
- [x] Electrostatic Matrix Slicing
- [ ] Ion Source Extraction Modeling
- [ ] 500 kV Accel Column Tracker
- [x] UI / Front-End Integration