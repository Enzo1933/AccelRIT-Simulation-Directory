#![allow(non_snake_case)]
use anyhow::Result;
use ndarray::{Array1, Array2, array};

use crate::{C_TM, MU0, PROTON_MASS};

/// Calculates the beam rigidity (B_rho)
/// Dimensions: T*m
/// Parameters: ke_mev [the kinetic energy in MeV]
fn beam_rigidity(ke_mev: f64) -> f64 {
    let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

    p / C_TM
}

/// Calculates the field gradient
/// Dimensions: T/m
/// Parameters: i [current], n [turns], r [radius]
fn field_gradient(i: f64, n: usize, r: f64, mu_r: f64) -> f64 {
    let ni = (n as f64) * i;
    let kappa = 1.0 / mu_r;

    (2.0 * MU0 * ni) / (r.powi(2) * (1.0 + kappa))
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

/// The envelope tracker struct
pub struct Tracker {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    x_f: f64, // final |x| position
    y_f: f64, // final |y| position
    total_length: f64,
    q1_end: f64,
    q2_start: f64,
    q2_end: f64,
    q3_start: f64,
    q3_end: f64,
    x_xover: Vec<f64>,
    y_xover: Vec<f64>,
    max_env_x: f64,
    max_env_y: f64,
}

impl Tracker {
    /// Track beam envelope through FDF triplet.
    /// Returns a Tracker data structure with z positions, x/y envelopes, region boundaries, crossovers, etc.
    pub fn new(
        L_mag_m: f64, // Magnet length in meters
        gap_m: f64,   // Gap length in meters
        drift_m: f64, // Drift length in meters
        g1: f64,      // The first field gradient
        g2: f64,      // The second field gradient
        energy_MeV: f64,
        x0: f64,
        xp0: f64,
        n_steps: usize,
    ) -> Result<Tracker> {
        let Brho = beam_rigidity(energy_MeV);
        let total_length = (3.0 * L_mag_m) + gap_m + (2.0 * drift_m);

        let regions = [
            ("quad", g1, L_mag_m),
            ("drift", 0.0, gap_m),
            ("quad", -g2, L_mag_m),
            ("drift", 0.0, drift_m),
            ("quad", -g1, L_mag_m),
            ("drift", 0.0, drift_m),
        ];

        let mut x = vec![x0];
        let mut y = vec![x0];
        let mut z = vec![0.0];

        let mut x_state: Array1<f64> = array![x0, xp0];
        let mut y_state: Array1<f64> = array![x0, xp0];

        for (r, g, length) in regions {
            let n = usize::max((n_steps as f64 * length / total_length) as usize, 4);
            let dz = length / n as f64;

            let (Mx, My) = match r {
                "quad" => quad_transfer_matrix(g, dz, Brho),
                _ => (drift_matrix(dz), drift_matrix(dz)),
            };

            for _ in 0..n {
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
        let q3_start = q2_end + drift_m;
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

    /// Performs Coordinate Descent to optimize both g1 and g2.
    pub fn optimize_triplet(
        L_mag_m: f64,
        gap_m: f64,
        drift_m: f64,
        energy_MeV: f64,
        x0: f64,
        xp0: f64,
        bore: f64,
    ) -> Option<(f64, f64)> {
        let mut g1 = 20.0; // Initial guess for outer quads
        let mut g2 = 20.0; // Initial guess for inner quad

        let mut last_g1 = g1;
        let mut last_g2 = g2;

        // Coordinate Descent: Toggle between optimizing g1 and g2
        for _ in 0..15 {
            // Optimize g1 (holding g2 constant)
            g1 = Self::golden_section_1d(0.1, 100.0, |test_g1| {
                Self::eval_fitness(
                    L_mag_m, gap_m, drift_m, test_g1, g2, energy_MeV, x0, xp0, bore,
                )
            });

            // Optimize g2 (holding g1 constant)
            g2 = Self::golden_section_1d(0.1, 100.0, |test_g2| {
                Self::eval_fitness(
                    L_mag_m, gap_m, drift_m, g1, test_g2, energy_MeV, x0, xp0, bore,
                )
            });

            last_g1 = g1;
            last_g2 = g2;
        }

        // Final validation
        let final_score =
            Self::eval_fitness(L_mag_m, gap_m, drift_m, g1, g2, energy_MeV, x0, xp0, bore);
        if final_score > 1.0 {
            None
        } else {
            Some((g1, g2))
        }
    }

    /// Generic 1D Golden Section Search
    fn golden_section_1d<F>(mut a: f64, mut b: f64, mut cost_fn: F) -> f64
    where
        F: FnMut(f64) -> f64,
    {
        let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
        let mut c = b - phi * (b - a);
        let mut d = a + phi * (b - a);

        for _ in 0..40 {
            if cost_fn(c) < cost_fn(d) {
                b = d;
            } else {
                a = c;
            }
            c = b - phi * (b - a);
            d = a + phi * (b - a);
        }
        (a + b) / 2.0
    }

    /// Fitness function 
    fn eval_fitness(
        L: f64,
        gap: f64,
        drift: f64,
        g1: f64,
        g2: f64,
        en: f64,
        x0: f64,
        xp0: f64,
        bore: f64,
    ) -> f64 {
        match Self::new(L, gap, drift, g1, g2, en, x0, xp0, 150) {
            Ok(t) => {
                if t.max_env_x > bore || t.max_env_y > bore {
                    return 1e9; // Penalty for hitting the pipe
                }
                // Balance symmetry and minimize final divergence/size
                let asymmetry = (t.x_f - t.y_f).abs();
                let size = t.x_f + t.y_f;
                (asymmetry / (size + 1e-10)) + size
            }
            Err(_) => 1e12,
        }
    }
}
