#![allow(non_snake_case)]
use anyhow::Result;
use ndarray::{Array2, array};

use crate::{C_Tm, PROTON_MASS};

/// Calculates the beam rigidity (B_rho)
/// Dimensions: T*m
/// Parameters: ke_mev [the kinetic energy in MeV]
pub fn beam_rigidity(ke_mev: f32) -> f32 {
    let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

    p / C_Tm
}

/// Calculates the quadrupole transfer matrix
pub fn quad_transfer_matrix(
    g: f32,     // Field gradient
    L: f32,     // Effective length
    B_rho: f32, // The beam rigidity
) -> Result<(Array2<f32>, Array2<f32>)> {
    let k2 = g / B_rho; // The magnetic gradient
    let k = k2.abs().sqrt();

    if k2.abs() < 1e-9 {
        // Near-zero gradient: both planes are drifts
        let drift = drift_matrix(L);
        return Ok((drift.clone(), drift));
    }

    if k2 > 0.0 {
        // Focusing in x, Defocusing in y
        let M_f = array![
            [(L * k.sqrt()).cos(), (L * k.sqrt()).sin() / k.sqrt()],
            [-1.0 * (L * k.sqrt()).sin() * k.sqrt(), (L * k.sqrt()).cos()]
        ];
        let M_d = array![
            [(L * k.sqrt()).cosh(), (L * k.sqrt()).sinh() / k.sqrt()],
            [
                -1.0 * (L * k.sqrt()).sinh() * k.sqrt(),
                (L * k.sqrt()).cosh()
            ]
        ];

        Ok((M_f, M_d))
    } else {
        // Defocusing in x, Focusing in y
        let M_f = array![
            [(L * k.sqrt()).cosh(), (L * k.sqrt()).sinh() / k.sqrt()],
            [
                -1.0 * (L * k.sqrt()).sinh() * k.sqrt(),
                (L * k.sqrt()).cosh()
            ]
        ];
        let M_d = array![
            [(L * k.sqrt()).cos(), (L * k.sqrt()).sin() / k.sqrt()],
            [-1.0 * (L * k.sqrt()).sin() * k.sqrt(), (L * k.sqrt()).cos()]
        ];

        Ok((M_f, M_d))
    }
}

/// Returns a drift matrix
pub fn drift_matrix(L: f32) -> Array2<f32> {
    array![[1.0, L], [0.0, 1.0]]
}
