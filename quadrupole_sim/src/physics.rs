#![allow(non_snake_case)]
use anyhow::Result;
use ndarray::{Array2, array};

use crate::{C_TM, MU0, PROTON_MASS};

/// Calculates the beam rigidity (B_rho)
/// Dimensions: T*m
/// Parameters: ke_mev [the kinetic energy in MeV]
pub fn beam_rigidity(ke_mev: f32) -> f32 {
    let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

    p / C_TM
}

/// Calculates the field gradient
/// Dimensions: T/m
/// Parameters: i [current], n [turns], r [radius]
pub fn field_gradient(i: f32, n: usize, r: f32, mu_r: f32) -> f32 {
    let ni = (n as f32) * i;
    let kappa = 1.0 / mu_r;  

    (2.0 * MU0 * ni) / (r.powi(2) * (1.0 + kappa))
}

/// Calculates the quadrupole transfer matrix
pub fn quad_transfer_matrix(
    g: f32,     // Field gradient
    L: f32,     // Effective length
    B_rho: f32, // The beam rigidity
) -> (Array2<f32>, Array2<f32>) {
    let k = g / B_rho; // Focusing strength
    let kr = k.abs().sqrt();

    if k.abs() < 1e-9 {
        // Near-zero gradient: both planes are drifts
        let drift = drift_matrix(L);
        return (drift.clone(), drift);
    }

    if k > 0.0 {
        // Focusing in x, Defocusing in y
        let M_f = array![
            [(L * kr).cos(), (L * kr).sin() / kr],
            [-1.0 * (L * kr).sin() * k.sqrt(), (L * kr).cos()]
        ];
        let M_d = array![
            [(L * kr).cosh(), (L * kr).sinh() / kr],
            [(L * kr).sinh() * k, (L * kr).cosh()]
        ];

        (M_f, M_d)
    } else {
        // Defocusing in x, Focusing in y
        let M_f = array![
            [(L * kr).cosh(), (L * kr).sinh() / kr],
            [-1.0 * (L * kr).sinh() * kr, (L * kr).cosh()]
        ];
        let M_d = array![
            [(L * kr).cos(), (L * kr).sin() / kr],
            [-1.0 * (L * kr).sin() * kr, (L * kr).cos()]
        ];

        (M_f, M_d)
    }
}

/// Returns a drift matrix
pub fn drift_matrix(L: f32) -> Array2<f32> {
    array![[1.0, L], [0.0, 1.0]]
}

/// The Envelope struct, represents a beam envelope
struct Envelope {
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    // ...
}

impl Envelope {
    /// Track beam envelope through FD doublet.
    /// Returns an Envelope data structure with z positions, x/y envelopes, region boundaries, crossovers, etc.
    pub fn track_envelope(
        bore_m: f32,  // Bore radius in meters
        L_mag_m: f32, // Magnet length in meters
        gap_m: f32,   // Gap length in meters
        drift_m: f32, // Drift length in meters
        g1: f32,      // The field gradient
        energy_MeV: f32,
        x0: f32,
        xp0: f32,
        mut n_steps: Option<usize>,
    ) -> Result<Envelope> {
        if n_steps == None {
            n_steps = Some(400); // Default steps is 400
        }

        let Brho = beam_rigidity(energy_MeV);
        let total_length = L_mag_m + gap_m + L_mag_m + drift_m;

        let regions = [
            ("quad", g1, L_mag_m),
            ("drift", 0.0, gap_m),
            ("quad", -g1, L_mag_m),
            ("drift", 0.0, drift_m),
        ];

        todo!()
    }
}
