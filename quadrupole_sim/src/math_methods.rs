use nalgebra::{SMatrix, Vector2, matrix, vector};

use crate::{
    beam_and_tracker::{Beam, Tracker},
    magnet::MagnetGeometry,
};

pub fn x_prime(state: Vector2<f64>, k: f64) -> Vector2<f64> {
    vector![
        state[1],      // x'
        -k * state[0], // x''
    ]
}

pub fn y_prime(state: Vector2<f64>, k: f64) -> Vector2<f64> {
    vector![
        state[1],     // y'
        k * state[0], // y''
    ]
}

pub fn rk4_step<F>(state: Vector2<f64>, z: f64, dz: f64, f: F) -> Vector2<f64>
where
    F: Fn(f64, Vector2<f64>) -> Vector2<f64>,
{
    let k1 = f(z, state);
    let k2 = f(z + 0.5 * dz, state + 0.5 * dz * k1);
    let k3 = f(z + 0.5 * dz, state + 0.5 * dz * k2);
    let k4 = f(z + dz, state + dz * k3);

    state + (dz / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

pub fn get_residuals_from_mmf(
    mmf1: f64,
    mmf2: f64,
    beam: &Beam,
    geo: &MagnetGeometry,
) -> Vector2<f64> {
    let g1 = geo.field_gradient(mmf1);
    let g2 = geo.field_gradient(mmf2);

    // Track the beam using 75 steps (as in your original code)
    let t = Tracker::new(beam, geo, g1, g2, 75).unwrap();

    // Residual 0: The final X position (Aiming for 0.0)
    // Residual 1: The final Y position (Aiming for 0.0)
    // No absolute values. No averages. Just pure, differentiable physics.
    vector![t.x_f, t.y_f]
}

pub fn jacobian(mmf1: f64, mmf2: f64, beam: &Beam, geo: &MagnetGeometry) -> SMatrix<f64, 2, 2> {
    // Dynamically scale epsilon to be exactly 0.1% of the current MMF.
    // The .max(10.0) ensures that if MMF gets very close to 0, we don't divide by zero.
    let eps1 = (mmf1 * 0.001).abs().max(10.0);
    let eps2 = (mmf2 * 0.001).abs().max(10.0);

    let r = get_residuals_from_mmf(mmf1, mmf2, beam, geo);

    // Perturb MMF1
    let r1 = get_residuals_from_mmf(mmf1 + eps1, mmf2, beam, geo);

    // Perturb MMF2
    let r2 = get_residuals_from_mmf(mmf1, mmf2 + eps2, beam, geo);

    // Calculate slopes (Partial Derivatives)
    matrix![
        (r1[0] - r[0]) / eps1, (r2[0] - r[0]) / eps2;
        (r1[1] - r[1]) / eps1, (r2[1] - r[1]) / eps2;
    ]
}
