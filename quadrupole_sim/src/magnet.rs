use std::f64::consts::PI;

use crate::MU0;

pub struct MagnetGeometry {
    pub bore: f64,             // Radius of the beam pipe (m)
    pub l_mag: f64,            // Physical length of the magnet (m)
    pub w_pole: f64,           // Width of the pole tip face (m)
    pub l_iron: f64,           // Average path length through the iron yoke (m)
    pub a_iron: f64,           // Average cross-sectional area of the iron (m^2)
    pub mu_i: f64,             // The initial permeability
    pub b_sat: f64,            // The magnetic saturation
    pub inter_magnet_gap: f64, // The gap between magnets
}

impl MagnetGeometry {
    /// Constructs a new magnet geometry
    pub fn new(
        bore: f64,
        l_mag: f64,
        w_pole: f64,
        l_iron: f64,
        a_iron: f64,
        mu_i: f64,
        b_sat: f64,
        inter_magnet_gap: f64,
    ) -> Self {
        Self {
            bore,
            l_mag,
            w_pole,
            l_iron,
            a_iron,
            mu_i,
            b_sat,
            inter_magnet_gap,
        }
    }

    /// Calculates the reluctances of the nodes in the magnet
    pub fn calculate_reluctances(&self, mu_eff: f64) -> (f64, f64, f64) {
        // Reluctance of the gap
        let R_gap = 2.0 * self.bore / (MU0 * self.l_mag * self.w_pole);

        let arg = 2.0 * self.w_pole / (2.0 * self.bore);
        let fringing_factor = 1.0 + (2.0 * self.bore / self.a_iron.sqrt()) * arg.ln();

        // Use the permeance of the leak to calculate the reluctance of the leaking flux
        let P_leak = MU0 * self.l_mag * (1.0 + self.w_pole / self.bore).ln() / PI;
        let R_leak = 1.0 / P_leak;

        // The reluctance of the iron
        let R_iron = self.l_iron / (MU0 * mu_eff * self.a_iron);

        (R_gap / fringing_factor, R_leak, R_iron)
    }

    /// Calculates the effective relative permeability using the Froelich-Kennelly model.
    fn frohlich_kennelly_mu(&self, b_field: f64) -> f64 {
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
        let mu_eff = self.frohlich_kennelly_mu(b_field);
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
        let a = self.a_iron * r_load / self.b_sat;
        let b_coeff = -(self.a_iron * (r_load + r_iron_linear) + mmf / self.b_sat);
        let c = mmf;

        // Solve using the Quadratic Formula
        let discriminant = b_coeff.powi(2) - 4.0 * a * c;

        if discriminant < 0.0 {
            panic!("Fatal: Negative discriminant when running solving for B pole")
        }

        let final_b_iron = ((-b_coeff - discriminant.sqrt()) / (2.0 * a)).clamp(0.0, self.b_sat);

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
        (2.0 * b) / self.bore
    }

    /// The enge multiplier function f(z)
    /// G(z) = G_0 * f(z)
    /// Params: z = Current z, z_0 = Effective edge
    pub fn enge_multiplier(&self, z: f64, z_0: f64) -> f64 {
        let a = [
            0.2965,  // a0: Constant offset
            2.2375,  // a1: Primary linear slope
            -0.0153, // a2: Curvature at the corner
            0.8847,  // a3: Higher order correction
            -0.2981, // a4: Tail smoothing
            0.0462,  // a5: Asymptotic damping
        ];

        let m = (z - z_0) / (2.0 * self.bore);

        1.0 / (1.0
            + (a[0]
                + a[1] * m
                + a[2] * m.powf(2.0)
                + a[3] * m.powf(3.0)
                + a[4] * m.powf(4.0)
                + a[5] * m.powf(5.0))
            .exp())
    }

    /// Returns raw gradient at position z accounting for fringe fields
    pub fn effective_gradient(
        &self,
        g0: f64,
        z: f64,
        z_entry: f64,
        z_exit: f64,
        L_eff: f64,
    ) -> f64 {
        let f_entry = self.enge_multiplier(z, z_entry);
        let f_exit = self.enge_multiplier(-(z - z_exit), 0.0);

        let f = f_entry * f_exit;
        let scale = self.l_mag / L_eff;

        g0 * f * scale
    }

    /// Effective magnetic length — integral of f(z) dz
    /// Numerically integrate the Enge function over a wide range
    pub fn effective_length(&self, z_entry: f64, z_exit: f64) -> f64 {
        let n = 500;
        let z_min = z_entry - 5.0 * self.bore;
        let z_max = z_exit + 5.0 * self.bore;
        let dz = (z_max - z_min) / n as f64;

        (0..n)
            .map(|i| {
                let z = z_min + i as f64 * dz;

                let f_entry = self.enge_multiplier(z, z_entry);
                let f_exit = self.enge_multiplier(-(z - z_exit), 0.0);

                (f_entry * f_exit) * dz
            })
            .sum()
    }
}
