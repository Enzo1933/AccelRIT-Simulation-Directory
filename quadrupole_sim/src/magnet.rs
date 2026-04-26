use std::f64::consts::PI;

use crate::MU0;

pub struct MagnetGeometry {
    pub r_gap: f64,  // Radius of the beam pipe (m)
    pub l_mag: f64,  // Physical length of the magnet (m)
    pub w_pole: f64, // Width of the pole tip face (m)
    pub l_iron: f64, // Average path length through the iron yoke (m)
    pub a_iron: f64, // Average cross-sectional area of the iron (m^2)
    pub mu_i: f64,   // The initial permeability
    pub b_sat: f64,  // The magnetic saturation
}

impl MagnetGeometry {
    /// Constructs a new magnet geometry
    pub fn new(
        r_gap: f64,
        l_mag: f64,
        w_pole: f64,
        l_iron: f64,
        a_iron: f64,
        mu_i: f64,
        b_sat: f64,
    ) -> Self {
        Self {
            r_gap,
            l_mag,
            w_pole,
            l_iron,
            a_iron,
            mu_i,
            b_sat,
        }
    }

    /// Calculates the reluctances of the nodes in the magnet
    fn calculate_reluctances(&self, mu_eff: f64) -> (f64, f64, f64) {
        // Reluctance of the gap
        let R_gap = self.r_gap / (MU0 * self.l_mag * self.w_pole);

        // Use the permeance of the leak to calculate the reluctance of the leaking flux
        let P_leak = MU0 * self.l_mag * (1.0 + self.w_pole / self.r_gap).ln() / PI;
        let R_leak = 1.0 / P_leak;

        // The reluctance of the iron
        let R_iron = self.l_iron / (MU0 * mu_eff * self.a_iron);

        (R_gap, R_leak, R_iron)
    }

    /// Calculates the effective relative permeability using the Froelich-Kennelly model.
    fn froelich_kennelly_mu(&self, b_field: f64) -> f64 {
        let b_abs = b_field.abs();

        if b_abs >= self.b_sat {
            return 1.0; // Fully saturated iron behaves like vacuum
        }

        let mu_eff = self.mu_i * (1.0 - (b_abs / self.b_sat));
        mu_eff.max(1.0)
    }

    /// Solves for the fluxes
    pub fn solve_fluxes(&self, b_field: f64, mmf: f64) -> (f64, f64, f64) {
        // Find the effective mu and the reluctances
        let mu_eff = self.froelich_kennelly_mu(b_field);
        let (r_gap, r_leak, r_iron) = self.calculate_reluctances(mu_eff);

        // Convert to permeance to get the load reluctance
        let p_load = 1.0 / r_gap + 1.0 / r_leak;
        let r_load = 1.0 / p_load;

        // Find total reluctance
        let R = r_load + r_iron;

        // Calculate phi
        let phi = mmf / R;
        let phi_gap = phi * (r_load / r_gap);
        let phi_leak = phi * (r_load / r_leak);

        (phi, phi_gap, phi_leak)
    }

    /// Magnetic Potential Energies
    pub fn magnetic_energies(&self, b_field: f64, mmf: f64) -> (f64, f64, f64) {
        let (phi, phi_gap, phi_leak) = self.solve_fluxes(b_field, mmf);
        let E = |phi: f64| 0.5 * phi * mmf;

        (E(phi), E(phi_gap), E(phi_leak))
    }

    /// Solves for the B field at the tip of a quadrupole
    pub fn solve_b_pole(&self, mmf: f64) -> f64 {
        // Get the baseline linear reluctances
        let (r_gap, r_leak, r_iron_linear) = self.calculate_reluctances(self.mu_i);
        let r_load = (r_gap * r_leak) / (r_gap + r_leak);

        // Setup the Quadratic Equation coefficients (a*B^2 + b*B + c = 0)
        let a = r_load / self.b_sat;
        let b_coeff = -(r_load + r_iron_linear + (mmf / (self.a_iron * self.b_sat)));
        let c = mmf / self.a_iron;

        // Solve using the Quadratic Formula
        let discriminant = b_coeff.powi(2) - 4.0 * a * c;

        if discriminant < 0.0 {
            return self.b_sat;
        }

        let final_b_iron = (-b_coeff + discriminant.sqrt()) / (2.0 * a);

        // Flux Divider Rule
        let phi_total = final_b_iron * self.a_iron;
        let phi_gap = phi_total * (r_load / r_gap);

        // Return the localized B-field at the pole tip
        phi_gap / (self.l_mag * self.w_pole)
    }

    /// Calculates the field gradient
    /// Dimensions: T/m
    pub fn field_gradient(&self, mmf: f64) -> f64 {
        let b = self.solve_b_pole(mmf);
        (2.0 * b) / self.r_gap
    }
}
