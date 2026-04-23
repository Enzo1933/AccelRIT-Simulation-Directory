#![allow(non_snake_case)]

use anyhow::{Ok, Result};
use ndarray::{Array1, Array2, array};
use std::fs::File;
use std::io::Write;

use crate::{C_TM, MU0, PROTON_MASS};

/// Calculates the beam rigidity (B_rho)
/// Dimensions: T*m
/// Parameters: ke_mev [the kinetic energy in MeV]
pub fn beam_rigidity(ke_mev: f64) -> f64 {
    let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

    p / C_TM
}

/// Solves for the B field of a quadrupole
/// Solve the ODE: F(B) = l_mat*H(B)+l_gap*(B/mu_0) - NI
pub fn solve_b_pole(i: f64, n: usize, r: f64, mu_r: f64, sat: f64, l_mag_m: f64, gap_m: f64) -> f64 {
    let ni = i * (n as f64);
    let mut b = (MU0 * ni) / r; // initial guess
    let epsilon = 1e-6;
    let mut run = true;

    while run {
        let old_b = b;

        let mu = effective_permeability(mu_r, sat, b);
        let b_target = (MU0 * mu * ni) / r;

        b = b * 0.7 + 0.3 * b_target;

        if (b - old_b).abs() < epsilon {
            run = false;
        }
    }

    todo!("Add the l_gap and l_mat parameters, run linear interpolation")
}

/// Calculates the field gradient
/// Dimensions: T/m
/// Parameters: i [current], n [turns], r [radius], mu_r [the relative permeability], sat [the saturation]
pub fn field_gradient(i: f64, n: usize, r: f64, mu_r: f64, sat: f64, l_mag: f64, l_gap: f64) -> f64 {
    let b = solve_b_pole(i, n, r, mu_r, sat, l_mag, l_gap);
    (2.0 * b) / r
}

/// Calculate Effective Permeability
fn effective_permeability(mu_r: f64, sat: f64, b_pole: f64) -> f64 {
    1.0 + (mu_r - 1.0) / (1.0 + (b_pole / sat).powi(4))
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
    pub L_mag_m: f64,    // Magnet length
    pub gap_m: f64,      // Inter-magnet gap (gap between quad-poles)
    pub drift_m: f64,    // Drift to the target
    pub energy_MeV: f64, // Kinetic energy
    pub x0: f64,         // x
    pub xp0: f64,        // x prime
}

impl Beam {
    pub fn new(L_mag_m: f64, gap_m: f64, drift_m: f64, energy_MeV: f64, x0: f64, xp0: f64) -> Self {
        Beam {
            L_mag_m,
            gap_m,
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
    pub x_f: f64, // final |x| position
    pub y_f: f64, // final |y| position
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
        g1: f64, // The first field gradient
        g2: f64, // The second field gradient
        n_steps: usize,
    ) -> Result<Tracker> {
        let L_mag_m: f64 = beam.L_mag_m;
        let gap_m: f64 = beam.gap_m;
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
            ("drift", 0.0, gap_m),   // Gap 2 (Was drift_m)
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

    /// Optimization using Newton-Raphson
    pub fn optimize_nr(
        args: &Beam,
        n1: usize,
        n2: usize,
        r: f64,
        mu_r: f64,
        sat: f64,
        l_mag_m: f64,
        gap_m: f64
    ) -> Option<(f64, f64)> {
        let mut i = array![20.0, 20.0];
        let eps = 1e-3; // Step size in Amps
        let learning_rate = 0.50;

        for _ in 0..50 {
            let res = Self::get_residuals_from_current(i[0], i[1], n1, n2, r, mu_r, sat, args, l_mag_m, gap_m);

            let res_i1 =
                Self::get_residuals_from_current(i[0] + eps, i[1], n1, n2, r, mu_r, sat, args, l_mag_m, gap_m);
            let res_i2 =
                Self::get_residuals_from_current(i[0], i[1] + eps, n1, n2, r, mu_r, sat, args, l_mag_m, gap_m);

            let jacobian = array![
                [(res_i1[0] - res[0]) / eps, (res_i2[0] - res[0]) / eps],
                [(res_i1[1] - res[1]) / eps, (res_i2[1] - res[1]) / eps]
            ];

            let det = jacobian[[0, 0]] * jacobian[[1, 1]] - jacobian[[0, 1]] * jacobian[[1, 0]];
            if det.abs() < 1e-14 {
                break;
            }

            let inv_j = array![
                [jacobian[[1, 1]] / det, -jacobian[[0, 1]] / det],
                [-jacobian[[1, 0]] / det, jacobian[[0, 0]] / det]
            ];

            let delta = &inv_j.dot(&(-1.0 * res));
            i += &(delta * learning_rate);

            if delta.dot(delta).sqrt() < 1e-6 {
                return Some((i[0], i[1]));
            }
        }
        Some((i[0], i[1]))
    }

    fn get_residuals_from_current(
        i1: f64,
        i2: f64,
        n1: usize,
        n2: usize,
        r: f64,
        mu_r: f64,
        sat: f64,
        beam: &Beam,
        l_mag_m: f64,
        gap_m: f64,
    ) -> Array1<f64> {
        let g1 = field_gradient(i1, n1, r, mu_r, sat, l_mag_m, gap_m);
        let g2 = field_gradient(i2, n2, r, mu_r, sat, l_mag_m, gap_m);

        // Check if we are saturating (B = G * r / 2)
        let b_pole2 = (g2 * r) / 2.0;
        let saturation_penalty = if b_pole2 > 1.8 {
            (b_pole2 - 1.8) * 1000.0
        } else {
            0.0
        };

        let t = Self::new(beam, g1, g2, 50).unwrap();

        // The optimizer now has to balance focusing with NOT saturating the iron
        array![t.x_f - t.y_f, t.x_f + saturation_penalty]
    }

    /// Exports the optimized profile as a CSV for IBSimu import.
    pub fn export_to_ibsimu(
        beam: &Beam,
        n1: usize,
        n2: usize,
        r: f64,
        mu_r: f64,
        sat: f64,
        l_mag_m: f64,
        gap_m: f64,
    ) -> Result<()> {
        let mut file = File::create("../beam_tracing.csv")?;
        let (i1, i2) = Self::optimize_nr(beam, n1, n2, r, mu_r, sat, l_mag_m, gap_m).unwrap();

        let g1 = field_gradient(i1, n1, r, mu_r, sat, l_mag_m, gap_m);
        let g2 = field_gradient(i2, n2, r, mu_r, sat, l_mag_m, gap_m);
        let final_tracker = Tracker::new(beam, g1, g2, 500)?;

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
    pub fn export_femm_lookup(
        beam: &Beam,
        n1: usize,
        n2: usize,
        r: f64,
        mu_r: f64,
        sat: f64,
        l_mag_m: f64,
        gap_m: f64,
    ) -> Result<()> {
        let (i1, i2) = Self::optimize_nr(beam, n1, n2, r, mu_r, sat, l_mag_m, gap_m).unwrap();

        let g1 = field_gradient(i1, n1, r, mu_r, sat, l_mag_m, gap_m);
        let g2 = field_gradient(i2, n2, r, mu_r, sat, l_mag_m, gap_m);
        let mut file = File::create("../FEMM-Lookup.csv")?;

        let b_pole1 = solve_b_pole(i1, n1, r, mu_r, sat, l_mag_m, gap_m);
        let mu_eff1 = effective_permeability(mu_r, sat, b_pole1);

        let b_pole2 = solve_b_pole(i2, n2, r, mu_r, sat, l_mag_m, gap_m);
        let mu_eff2 = effective_permeability(mu_r, sat, b_pole2);

        writeln!(
            file,
            "Magnet,Gradient(T/m),Current(A),Turns,Mu_eff,Radius(m)\n\
             Outer_Quads,{},{},{},{},{}\n\
             Inner_Quad,{},{},{},{},{}",
            g1, i1, n1, mu_eff1, r, g2, i2, n2, mu_eff2, r
        )?;

        todo!()
    }
}
