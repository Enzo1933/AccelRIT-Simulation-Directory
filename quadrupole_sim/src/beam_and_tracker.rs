#![allow(non_snake_case)]

use anyhow::{Ok, Result};
use ndarray::{Array1, Array2, array};
use std::fs::File;
use std::io::Write;

use crate::{C_TM, PROTON_MASS, magnet::MagnetGeometry};

/// Calculates the beam rigidity (B_rho)
/// Dimensions: T*m
/// Parameters: ke_mev [the kinetic energy in MeV]
pub fn beam_rigidity(ke_mev: f64) -> f64 {
    let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

    p / C_TM
}

/// Calculates the quadrupole transfer matrix
fn quad_transfer_matrix(
    g: f64,     // Field gradient
    L: f64,     // Effective length
    B_rho: f64, // The beam rigidity
) -> (Array2<f64>, Array2<f64>) {
    let k = g / B_rho; // Focusing strength
    let kr = k.abs().sqrt();

    if k.abs() < 1e-9 {
        // Near-zero gradient: both planes are drifts
        let drift = drift_matrix(L);
        return (drift.clone(), drift);
    }

    if k > 0.0 {
        // Focusing in x, Defocusing in y
        let M_x = array![
            [(L * kr).cos(), (L * kr).sin() / kr],
            [-1.0 * (L * kr).sin() * kr, (L * kr).cos()]
        ];
        let M_y = array![
            [(L * kr).cosh(), (L * kr).sinh() / kr],
            [(L * kr).sinh() * kr, (L * kr).cosh()]
        ];

        (M_x, M_y)
    } else {
        // Defocusing in x, Focusing in y
        let M_x = array![
            [(L * kr).cosh(), (L * kr).sinh() / kr],
            [1.0 * (L * kr).sinh() * kr, (L * kr).cosh()]
        ];
        let M_y = array![
            [(L * kr).cos(), (L * kr).sin() / kr],
            [-1.0 * (L * kr).sin() * kr, (L * kr).cos()]
        ];

        (M_x, M_y)
    }
}

/// Returns a drift matrix
fn drift_matrix(L: f64) -> Array2<f64> {
    array![[1.0, L], [0.0, 1.0]]
}

fn find_crossovers(arr: &[f64], z: &[f64]) -> Vec<f64> {
    let mut crossovers = Vec::new();
    for i in 1..arr.len() {
        if arr[i - 1] * arr[i] < 0.0 {
            let frac = arr[i - 1].abs() / (arr[i - 1].abs() + arr[i].abs());
            crossovers.push(z[i - 1] + frac * (z[i] - z[i - 1]));
        }
    }
    crossovers
}

/// Beam struct
pub struct Beam {
    pub drift_m: f64,    // Drift to the target
    pub energy_MeV: f64, // Kinetic energy
    pub x0: f64,         // x
    pub xp0: f64,        // x prime
}

impl Beam {
    pub fn new(drift_m: f64, energy_MeV: f64, x0: f64, xp0: f64) -> Self {
        Beam {
            drift_m,
            energy_MeV,
            x0,
            xp0,
        }
    }
}

/// The envelope tracker struct
pub struct Tracker {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub x_f: f64,
    pub y_f: f64,
    pub total_length: f64,
    pub q1_end: f64,
    pub q2_start: f64,
    pub q2_end: f64,
    pub q3_start: f64,
    pub q3_end: f64,
    pub x_xover: Vec<f64>,
    pub y_xover: Vec<f64>,
    pub max_env_x: f64,
    pub max_env_y: f64,
}

impl Tracker {
    /// Track beam envelope through FDF triplet.
    /// Returns a Tracker data structure with z positions, x/y envelopes, region boundaries, crossovers, etc.
    pub fn new(
        beam: &Beam,
        geo: &MagnetGeometry,
        g1: f64,
        g2: f64,
        n_steps: usize,
    ) -> Result<Tracker> {
        let L_mag_m: f64 = geo.l_mag;
        let gap_m: f64 = geo.r_gap;
        let drift_m: f64 = beam.drift_m;
        let energy_MeV = beam.energy_MeV;
        let x0 = beam.x0;
        let xp0 = beam.xp0;

        let Brho = beam_rigidity(energy_MeV);
        let total_length = (3.0 * L_mag_m) + (2.0 * gap_m) + drift_m;

        let regions = [
            ("quad", g1, L_mag_m),   // Q1
            ("drift", 0.0, gap_m),   // Gap 1
            ("quad", -g2, L_mag_m),  // Q2
            ("drift", 0.0, gap_m),   // Gap 2
            ("quad", g1, L_mag_m),   // Q3 (Must match Q1 polarity)
            ("drift", 0.0, drift_m), // Final Drift to target
        ];

        let mut x = vec![x0];
        let mut y = vec![x0];
        let mut z = vec![0.0];

        let mut x_state: Array1<f64> = array![x0, xp0];
        let mut y_state: Array1<f64> = array![x0, xp0];

        for (r, g, length) in regions {
            let n = usize::max((n_steps as f64 * length / total_length) as usize, 4);
            let dz = length / n as f64;

            for _ in 0..n {
                let (Mx, My) = match r {
                    "quad" => quad_transfer_matrix(g, dz, Brho),
                    _ => (drift_matrix(dz), drift_matrix(dz)),
                };

                x_state = Mx.dot(&x_state);
                y_state = My.dot(&y_state);

                z.push(z.last().unwrap() + dz);
                x.push(x_state[0]);
                y.push(y_state[0]);
            }
        }

        let x_f = x_state[0].abs();
        let y_f = y_state[0].abs();

        let max_env_x = x.iter().map(|v| v.abs()).fold(f64::NEG_INFINITY, f64::max);
        let max_env_y = y.iter().map(|v| v.abs()).fold(f64::NEG_INFINITY, f64::max);

        let x_xover = find_crossovers(&x, &z);
        let y_xover = find_crossovers(&y, &z);

        let q1_end = L_mag_m;
        let q2_start = q1_end + gap_m;
        let q2_end = q2_start + L_mag_m;
        let q3_start = q2_end + gap_m;
        let q3_end = q3_start + L_mag_m;

        Ok(Tracker {
            x,
            y,
            z,
            x_f,
            y_f,
            total_length,
            q1_end,
            q2_start,
            q2_end,
            q3_start,
            q3_end,
            x_xover,
            y_xover,
            max_env_x,
            max_env_y,
        })
    }

    /// Optimize magneto-motive force
    pub fn optimize_mmf(beam: &Beam, geo: &MagnetGeometry) -> Option<(f64, f64)> {
        let mut mmf1 = 15000.0; // Outer quads
        let mut mmf2 = 15000.0; // Inner quad

        let mut run = true;
        let eps = 10.0;
        let lambda = 0.5; // Damping factor to prevent overshoot

        while run {
            let (res_asym, res_size) = Tracker::get_residuals_from_mmf(mmf1, mmf2, beam, geo);

            // Convergence Check 
            if res_asym.abs() < 1e-6 && res_size.abs() < 1e-6 {
                run = false;
            }

            let (nudge1_asym, nudge1_size) =
                Tracker::get_residuals_from_mmf(mmf1 + eps, mmf2, beam, geo);
            let j11 = (nudge1_asym - res_asym) / eps; // d(Asymmetry) / d(MMF1)
            let j21 = (nudge1_size - res_size) / eps; // d(Spot Size) / d(MMF1)

            let (nudge2_asym, nudge2_size) =
                Tracker::get_residuals_from_mmf(mmf1, mmf2 + eps, beam, geo);
            let j12 = (nudge2_asym - res_asym) / eps; // d(Asymmetry) / d(MMF2)
            let j22 = (nudge2_size - res_size) / eps; // d(Spot Size) / d(MMF2)

            let det = (j11 * j22) - (j12 * j21);

            if det.abs() < 1e-14 {
                run = false
            }

            // Calculate the raw Newton steps
            let delta_mmf1 = (j22 * res_asym - j12 * res_size) / det;
            let delta_mmf2 = (-j21 * res_asym + j11 * res_size) / det;

            // Apply the steps with the damping factor
            mmf1 -= lambda * delta_mmf1;
            mmf2 -= lambda * delta_mmf2;
        }

        Some((mmf1, mmf2))
    }

    fn get_residuals_from_mmf(
        mmf1: f64,
        mmf2: f64,
        beam: &Beam,
        geo: &MagnetGeometry,
    ) -> (f64, f64) {
        let g1 = geo.field_gradient(mmf1);
        let g2 = geo.field_gradient(mmf2);

        let t = Self::new(beam, geo, g1, g2, 50).unwrap();

        (t.x_f - t.y_f, t.x_f)
    }

    /// Exports the optimized profile as a CSV for IBSimu import.
    pub fn export_to_ibsimu(beam: &Beam, geo: &MagnetGeometry) -> Result<()> {
        let mut file = File::create("../beam_tracing.csv")?;
        let (mmf1, mmf2) = Self::optimize_mmf(beam, geo).unwrap();

        let g1 = geo.field_gradient(mmf1);
        let g2 = geo.field_gradient(mmf2);
        let final_tracker = Tracker::new(beam, geo, g1, g2, 500)?;

        writeln!(file, "z,x_env,y_env")?;
        for i in 0..final_tracker.z.len() {
            writeln!(
                file,
                "{},{},{}",
                final_tracker.z[i],
                final_tracker.x[i].abs(),
                final_tracker.y[i].abs()
            )?;
        }
        Ok(())
    }

    /// Generates a CSV lookup table for FEMM import
    pub fn export_femm_lookup(beam: &Beam, geo: &MagnetGeometry) -> Result<()> {
        let (mmf1, mmf2) = Self::optimize_mmf(beam, geo).unwrap();

        let g1 = geo.field_gradient(mmf1);
        let g2 = geo.field_gradient(mmf2);
        let mut file = File::create("../FEMM-Lookup.csv")?;

        writeln!(
            file,
            "Magnet,Gradient(T/m),Magnetomotive Force (A*t)\n\
             Outer_Quads,{},{}\n\
             Inner_Quad,{},{}",
            g1, mmf1, g2, mmf2
        )?;

        Ok(())
    }
}
