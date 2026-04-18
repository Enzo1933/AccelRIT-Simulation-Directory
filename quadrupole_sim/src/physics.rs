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
pub fn field_gradient(i: f32, n: usize, r: f32) -> f32 {
    let ampere_turns = (n as f32) * i;

    2.0 * MU0 * ampere_turns / r.powi(2)
}

/// Calculates the quadrupole transfer matrix
pub fn quad_transfer_matrix(
    g: f32,     // Field gradient
    L: f32,     // Effective length
    B_rho: f32, // The beam rigidity
) -> (Array2<f32>, Array2<f32>) {
    let k2 = g / B_rho; // Magnetic field strength
    let k = k2.abs().sqrt();

    if k2.abs() < 1e-9 {
        // Near-zero gradient: both planes are drifts
        let drift = drift_matrix(L);
        return (drift.clone(), drift);
    }

    if k2 > 0.0 {
        // Focusing in x, Defocusing in y
        let M_f = array![
            [(L * k).cos(), (L * k).sin() / k],
            [-1.0 * (L * k).sin() * k.sqrt(), (L * k).cos()]
        ];
        let M_d = array![
            [(L * k).cosh(), (L * k).sinh() / k],
            [-1.0 * (L * k).sinh() * k, (L * k).cosh()]
        ];

        (M_f, M_d)
    } else {
        // Defocusing in x, Focusing in y
        let M_f = array![
            [(L * k).cosh(), (L * k).sinh() / k],
            [-1.0 * (L * k).sinh() * k, (L * k).cosh()]
        ];
        let M_d = array![
            [(L * k).cos(), (L * k).sin() / k],
            [-1.0 * (L * k).sin() * k, (L * k).cos()]
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
        bore_m: f32,
        L_mag_m: f32,
        gap_m: f32,
        drift_m: f32,
        g1: f32,
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
