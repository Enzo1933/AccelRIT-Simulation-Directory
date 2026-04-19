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
        let M_f = array![
            [(L * kr).cos(), (L * kr).sin() / kr],
            [-1.0 * (L * kr).sin() * kr, (L * kr).cos()]
        ];
        let M_d = array![
            [(L * kr).cosh(), (L * kr).sinh() / kr],
            [(L * kr).sinh() * kr, (L * kr).cosh()]
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
    x_xover: Vec<f64>,
    y_xover: Vec<f64>,
    max_env_x: f64,
    max_env_y: f64,
}

impl Tracker {
    /// Track beam envelope through FD doublet.
    /// Returns a Tracker data structure with z positions, x/y envelopes, region boundaries, crossovers, etc.
    pub fn new(
        L_mag_m: f64, // Magnet length in meters
        gap_m: f64,   // Gap length in meters
        drift_m: f64, // Drift length in meters
        g1: f64,      // The field gradient
        energy_MeV: f64,
        x0: f64,
        xp0: f64,
        n_steps: usize,
    ) -> Result<Tracker> {
        let Brho = beam_rigidity(energy_MeV);
        let total_length = L_mag_m + gap_m + L_mag_m + drift_m;

        let regions = [
            ("quad", g1, L_mag_m),
            ("drift", 0.0, gap_m),
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

            for _ in 0..n {
                match r {
                    "quad" => {
                        let (Mx, My) = quad_transfer_matrix(g, dz, Brho);

                        x_state = Mx.dot(&x_state);
                        y_state = My.dot(&y_state);
                    }
                    _ => {
                        x_state = drift_matrix(dz).dot(&x_state);
                        y_state = drift_matrix(dz).dot(&y_state);
                    }
                }

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

        Ok(Tracker {
            x,
            y,
            z,
            x_f,
            y_f,
            total_length,
            q1_end: L_mag_m,
            q2_start: L_mag_m + gap_m,
            q2_end: 2.0 * L_mag_m + gap_m,
            x_xover,
            y_xover,
            max_env_x,
            max_env_y,
        })
    }

    /// Uses Golden Section Search for gradient optimization
    pub fn optimize_gradient(&self, bore_m: f64) -> Option<f64> {
        let eval = |g: f64| -> Option<f64> {
            if self.max_env_x > bore_m * 0.95 {
                return None;
            }
            if self.max_env_y > bore_m * 0.95 {
                return None;
            }

            let avg = (self.x_f + self.y_f) / 2.0;
            if avg < 1e-12 {
                return None;
            }

            Some((self.x_f - self.y_f).abs() / avg)
        };

        let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
        let tol = 1e-4;

        let mut a = 0.1_f64;
        let mut b = 50.0_f64;
        let mut c = b - phi * (b - a);
        let mut d = a + phi * (b - a);

        while (b - a).abs() > tol {
            let fc = eval(c).unwrap_or(f64::INFINITY);
            let fd = eval(d).unwrap_or(f64::INFINITY);

            if fc < fd {
                b = d;
            } else {
                a = c;
            }

            c = b - phi * (b - a);
            d = a + phi * (b - a);
        }

        let g_opt = (a + b) / 2.0;

        // Verify the solution is actually valid before returning
        eval(g_opt).map(|_| g_opt)
    }
}
