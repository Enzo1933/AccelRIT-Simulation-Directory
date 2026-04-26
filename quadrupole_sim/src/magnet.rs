use std::f64::consts::{PI, SQRT_2};

use crate::MU0;

pub struct MagnetGeometry {
    pub r_gap: f64,  // Radius of the beam pipe (m)
    pub l_mag: f64,  // Physical length of the magnet (m)
    pub w_pole: f64, // Width of the pole tip face (m)
    pub h_pole: f64, // Height of the pole piece (m)
    pub l_iron: f64, // Average path length through the iron yoke (m)
    pub a_iron: f64, // Average cross-sectional area of the iron (m^2)
}

impl MagnetGeometry {
    /// Constructs a new magnet geometry
    pub fn new(
        r_gap: f64,  // Radius of the beam pipe (m)
        l_mag: f64,  // Physical length of the magnet (m)
        w_pole: f64, // Width of the pole tip face (m)
        h_pole: f64, // Height of the pole piece (m)
        l_iron: f64, // Average path length through the iron yoke (m)
        a_iron: f64, // Average cross-sectional area of the iron (m^2)
    ) -> Self {
        Self {
            r_gap,
            l_mag,
            w_pole,
            h_pole,
            l_iron,
            a_iron,
        }
    }

    /// Calculates the reluctances of the nodes in the magnet
    pub fn calculate_reluctances(&self, mu_eff: f64) -> (f64, f64, f64) {
        // Reluctance of the gap
        let R_gap = self.r_gap / (MU0 * self.l_mag * self.w_pole);

        // Use the permeance of the leak to calculate the reluctance of the leaking flux
        let P_leak = MU0 * self.l_mag * (1.0 + self.w_pole / self.r_gap).ln() / PI; 
        let R_leak = 1.0 / P_leak;

        // The reluctance of the iron
        let R_iron = self.l_iron / MU0 * mu_eff * self.a_iron;

        (R_gap, R_leak, R_iron)
    }
}

/// Solves for the B field of a quadrupole
pub fn solve_b_pole(i: f64, n: usize, l_iron: f64, l_gap: f64) -> f64 {
    todo!()
}

/// Calculates the field gradient
/// Dimensions: T/m
pub fn field_gradient(
    i: f64,   // Current
    n: usize, // Turns
    r: f64,   // Bore radius
    l_iron: f64,
    l_gap: f64,
) -> f64 {
    let b = solve_b_pole(i, n, l_iron, l_gap);
    (2.0 * b) / r
}
