#![allow(non_snake_case)]

use anyhow::{Ok, Result};
use nalgebra::{Vector2, vector};
use std::fs::File;
use std::io::Write;

use crate::{
    C_TM, PROTON_MASS,
    magnet::MagnetGeometry,
    math_methods::{get_residuals_from_mmf, jacobian, rk4_step, x_prime, y_prime},
};

/// Calculates the beam rigidity (B_rho)
/// Dimensions: T*m
/// Parameters: ke_mev [the kinetic energy in MeV]
pub fn beam_rigidity(ke_mev: f64) -> f64 {
    let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

    p / C_TM
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
        let gap_m: f64 = geo.inter_magnet_gap;
        let drift_m: f64 = beam.drift_m;
        let energy_MeV = beam.energy_MeV;
        let x0 = beam.x0;
        let xp0 = beam.xp0;

        let q1_start = 0.0;
        let q1_end = L_mag_m;

        let q2_start = L_mag_m + gap_m;
        let q2_end = L_mag_m + gap_m + L_mag_m;

        let q3_start = 2.0 * L_mag_m + 2.0 * gap_m;
        let q3_end = 3.0 * L_mag_m + 2.0 * gap_m;

        let Brho = beam_rigidity(energy_MeV);
        let total_length = (3.0 * L_mag_m) + (2.0 * gap_m) + drift_m;

        let L_eff_q1 = geo.effective_length(q1_start, q1_end);
        let L_eff_q2 = geo.effective_length(q2_start, q2_end);
        let L_eff_q3 = geo.effective_length(q3_start, q3_end);

        let regions = [
            ("quad", g1, L_mag_m, q1_start, q1_end, L_eff_q1),
            ("drift", 0.0, gap_m, 0.0, 0.0, 0.0),
            ("quad", -g2, L_mag_m, q2_start, q2_end, L_eff_q2),
            ("drift", 0.0, gap_m, 0.0, 0.0, 0.0),
            ("quad", g1, L_mag_m, q3_start, q3_end, L_eff_q3),
            ("drift", 0.0, drift_m, 0.0, 0.0, 0.0),
        ];

        let mut x = vec![x0];
        let mut y = vec![x0];
        let mut z = vec![0.0];

        let mut x_state = vector![x0, xp0];
        let mut y_state = vector![x0, xp0];

        for (r, g, length, z_entry, z_exit, L_eff) in regions {
            let n = usize::max((n_steps as f64 * length / total_length) as usize, 4);
            let dz = length / n as f64;

            for _ in 0..n {
                let z_curr = *z.last().unwrap();

                match r {
                    "quad" => {
                        let g_eff = geo.effective_gradient(g, z_curr, z_entry, z_exit, L_eff);
                        let k = g_eff / Brho;

                        x_state = rk4_step(x_state, z_curr, dz, |_z, s| x_prime(s, k));
                        y_state = rk4_step(y_state, z_curr, dz, |_z, s| y_prime(s, k));
                    }
                    _ => {
                        // Drift — k = 0, straight line
                        x_state = rk4_step(x_state, z_curr, dz, |_z, s| x_prime(s, 0.0));
                        y_state = rk4_step(y_state, z_curr, dz, |_z, s| y_prime(s, 0.0));
                    }
                }

                z.push(z_curr + dz);
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
    /// TODO: Fix this stuff
    pub fn optimize_mmf(beam: &Beam, geo: &MagnetGeometry) -> Option<(f64, f64)> {
        let mmf = Self::calculate_realistic_guess(beam, geo);
        let mut mmf1 = mmf[0];
        let mut mmf2 = mmf[1];

        // let lambda = 0.075; // damping
        let max_iter = 250;

        for i in 0..max_iter {
            let r = get_residuals_from_mmf(mmf1, mmf2, beam, geo);
            let r1 = r[0];
            let r2 = r[1];

            let cost = r1 * r1 + r2 * r2;
            println!("Iter {i:2} | MMF1: {mmf1:.1} MMF2: {mmf2:.1} | Cost: {cost:.3e}");

            if cost < 1e-8 {
                println!("Converged at iter {i}");
                break;
            }

            let j = jacobian(mmf1, mmf2, beam, geo);
            let det = j.determinant();
            if det.abs() < 1e-12 {
                break; // Safety net
            }

            let inv_j = j.try_inverse().unwrap();
            let mut delta = inv_j * r;

            // Never allow a step larger than 20% of the current MMF.
            // Stops the "teleportation" into saturation walls.
            let max_step1 = mmf1 * 0.20;
            let max_step2 = mmf2 * 0.20;

            delta[0] = delta[0].clamp(-max_step1, max_step1);
            delta[1] = delta[1].clamp(-max_step2, max_step2);

            // 2. APPLY THE UPDATE (with a modest lambda damping, e.g., 0.5)
            let lambda = 0.5;
            mmf1 = (mmf1 - delta[0] * lambda).clamp(10.0, 250_000.0);
            mmf2 = (mmf2 - delta[1] * lambda).clamp(10.0, 250_000.0);
        }

        let ratio = mmf2 / mmf1;

        if ratio > 2.0 || ratio < 1.5 {
            eprintln!("mmf2/mmf1 {ratio:?}, perhaps reconsider your options")
        }

        println!("Final: MMF1={mmf1:.1} MMF2={mmf2:.1}");
        Some((mmf1, mmf2))
    }

    /// Generates a physically realistic starting guess for MMF based on optical focal length
    pub fn calculate_realistic_guess(beam: &Beam, geo: &MagnetGeometry) -> Vector2<f64> {
        // Thin Lens Approximation
        let required_focal_length = beam.drift_m;

        // Optical strength (k) = 1 / (f * Effective Length)
        let k_estimate = 1.0 / (required_focal_length * geo.l_mag);

        // k = g / Brho  =>  g = k * Brho
        let g_estimate = k_estimate * beam_rigidity(beam.energy_MeV);
        // g = (2 * B_pole) / r  =>  B_pole = (g * r) / 2
        let b_pole_estimate = (g_estimate * geo.bore) / 2.0;

        // MMF = Flux * Reluctance
        let (r_gap, r_leak, r_iron) = geo.calculate_reluctances(geo.mu_i);
        let r_load = (r_gap * r_leak) / (r_gap + r_leak);
        let r_total = r_load + r_iron;

        let flux_total = b_pole_estimate * geo.a_iron;
        let base_mmf = flux_total * r_total;

        // Triplet Ratio Rule (Inner quad is ~1.6x stronger)
        let mmf1_guess = base_mmf;
        let mmf2_guess = base_mmf * 1.6;

        println!(
            "Physics Engine Estimated Starting MMFs: ({:.1}, {:.1})",
            mmf1_guess, mmf2_guess
        );

        vector![mmf1_guess, mmf2_guess]
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
