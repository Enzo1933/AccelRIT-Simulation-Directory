use crate::math_methods::sech;

pub struct EinzelGeometry {
    pub U_mid: f64, // Sitting voltage of the middle electrode
    pub L_mid: f64, // Length of the middle electrode
    pub R: f64,     // Radius of the three cylinders
}

impl EinzelGeometry {
    /// Constructs a new einzel lens
    pub fn new(U_mid: f64, L_mid: f64, R: f64) -> Self {
        Self { U_mid, L_mid, R }
    }

    /// Returns the potential at a certain z value
    pub fn voltage(&self, z: f64) -> f64 {
        let k = 1.318 / self.R;
        let u1 = k * (z + 0.5 * self.L_mid);
        let u2 = k * (z - 0.5 * self.L_mid);

        (self.U_mid / 2.0) * (u1.tanh() - u2.tanh())
    }

    /// Returns the E field at a certain z value
    pub fn e_field(&self, z: f64) -> f64 {
        let k = 1.318 / self.R;
        let u1 = k * (z + 0.5 * self.L_mid);
        let u2 = k * (z - 0.5 * self.L_mid);

        (self.U_mid / 2.0) * (sech(u1).powi(2) - sech(u2).powi(2))
    }

    /// Returns the radial focusing gradient E' at a certain z value
    pub fn rad_focusing_gradient(&self, z: f64) -> f64 {
        let k = 1.318 / self.R;
        let u1 = k * (z + 0.5 * self.L_mid);
        let u2 = k * (z - 0.5 * self.L_mid);

        -self.U_mid * k.powi(2) * (sech(u1).powi(2) * u1.tanh() - sech(u2).powi(2) * u2.tanh())
    }
}
