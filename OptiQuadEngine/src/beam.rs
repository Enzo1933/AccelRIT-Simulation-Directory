use crate::{C_TM, PROTON_MASS};

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

    /// Calculates the beam rigidity (B_rho)
    /// Dimensions: T*m
    /// Parameters: ke_mev [the kinetic energy in MeV]
    pub fn beam_rigidity(&self) -> f64 {
        let ke_mev = self.energy_MeV;
        let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

        p / C_TM
    }
}
