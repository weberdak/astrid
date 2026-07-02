use crate::paf::calc_paf_coords;
use crate::params::Params;
use crate::peak::{calc_dipolar_coupling, calc_shift};
use nalgebra::Matrix3;
use std::f64::consts::PI;

pub struct Residue {
    /// Index relative to reference residue in the sequence
    pub index: i16,
    /// Rotation angle of the residue in radians (0 - 2π)
    pub rotation: f64,
    /// Single-letter name of amino acid residue (e.g., "A", "C", "D", etc.)
    pub label: Option<[u8; 1]>,
    /// Simulated 15N chemical shift for the residue in ppm
    pub sim_shift: Option<f64>,
    /// Simulated 15N-1H dipolar coupling for the residue in kHz
    pub sim_dipolar_coupling: Option<f64>,
    /// Experimental 15N chemical shift for the residue in ppm
    pub exp_shift: Option<f64>,
    /// Experimental 15N-1H dipolar coupling for the residue in kHz
    pub exp_dipolar_coupling: Option<f64>,
    /// RMSD between simulated and experimental peaks for the residue
    pub rmsd: Option<f64>,
}

impl Residue {
    /// Create a new Residue instance with default values
    pub fn new(index: i16, params: &Params) -> Self {
        // 3.6 residues per turn in an ideal helix (100 degrees per residue)
        let rotation = (params.rotation - index as f64 * 100.0_f64.to_radians())
            .rem_euclid(2.0 * PI);
        Self {
            index,
            rotation,
            label: None,
            sim_shift: None,
            sim_dipolar_coupling: None,
            exp_shift: None,
            exp_dipolar_coupling: None,
            rmsd: None,
        }
    }

    /// Update the simulated shift and dipolar coupling for the residue
    pub fn calc_peak(&mut self, params: &Params, matrix_a: &Matrix3<f64>) {
        // 3.6 residues per turn in an ideal helix (100 degrees per residue)
        let paf_coords = calc_paf_coords(self.rotation, params, matrix_a);
        self.sim_shift = Some(calc_shift(params, &paf_coords));
        self.sim_dipolar_coupling =
            Some(calc_dipolar_coupling(params, &paf_coords).abs());
        self.calc_rmsd();
    }

    /// Update the experimental shift for the residue and recalculate RMSD
    pub fn set_exp_shift(&mut self, exp_shift: f64) {
        self.exp_shift = Some(exp_shift);
        self.calc_rmsd();
    }

    /// Unset the experimental shift for the residue and recalculate RMSD
    pub fn unset_exp_shift(&mut self) {
        self.exp_shift = None;
        self.calc_rmsd();
    }

    /// Update the experimental dipolar coupling for the residue and recalculate RMSD
    pub fn set_exp_dipolar_coupling(&mut self, exp_dipolar_coupling: f64) {
        self.exp_dipolar_coupling = Some(exp_dipolar_coupling);
        self.calc_rmsd();
    }

    /// Unset the experimental dipolar coupling for the residue and recalculate RMSD
    pub fn unset_exp_dipolar_coupling(&mut self) {
        self.exp_dipolar_coupling = None;
        self.calc_rmsd();
    }

    /// Update the label for the residue
    pub fn set_label(&mut self, label: [u8; 1]) {
        self.label = Some(label);
    }

    /// Print the residue information in a human-readable format
    pub fn print_info(&self) {
        println!(
            "Residue {}: rotation = {:.2} rad, sim_shift = {:?}, sim_dipolar_coupling = {:?}, exp_shift = {:?}, exp_dipolar_coupling = {:?}, rmsd = {:?}",
            self.index,
            self.rotation,
            self.sim_shift,
            self.sim_dipolar_coupling,
            self.exp_shift,
            self.exp_dipolar_coupling,
            self.rmsd
        );
    }

    /// Calculate the RMSD between simulated and experimental peaks for the residue
    fn calc_rmsd(&mut self) {
        // Scale coupling difference to match shift dispersion
        let scaling_factor = 10.0;

        let shift_diff_sq = self
            .sim_shift
            .zip(self.exp_shift)
            .map(|(sim, exp)| (sim - exp).powi(2));

        let coupling_diff_sq = self
            .sim_dipolar_coupling
            .zip(self.exp_dipolar_coupling)
            .map(|(sim, exp)| ((sim - exp) * scaling_factor).powi(2));

        self.rmsd = match (shift_diff_sq, coupling_diff_sq) {
            (Some(ssq), Some(csq)) => Some(((ssq + csq) / 2.0).sqrt()),
            (Some(ssq), None) => Some(ssq.sqrt()),
            (None, Some(csq)) => Some(csq.sqrt()),
            (None, None) => None,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix_a::calc_matrix_a;

    #[test]
    /// Test using SLN example from Weber et al. (2020)
    /// https://github.com/weberdak/pisa.py/tree/master/examples/sarcolipin/sln_explore_log.dat
    fn test_residue_calcs() {
        let mut params = Params::new();
        params.update_flip_angle(90.0);
        params.update_order_parameter(0.9);
        params.tilt = 24.6_f64.to_radians();
        params.rotation = 46.0_f64.to_radians();
        let mut residue = Residue::new(-3, &params);

        // RMSD must be empty before any calculations
        assert!(residue.rmsd.is_none());

        // Calculate simulated peak values for the residue, RMSD should be empty
        residue.calc_peak(&params, &calc_matrix_a(&params));
        assert!(residue.rmsd.is_none());

        // Sim chemical shift should be within 0.1 of 92.5
        assert!((residue.sim_shift.unwrap() - 92.52).abs() < 0.1);

        // Sim dipolar coupling should be within 0.01 of 4.52
        assert!((residue.sim_dipolar_coupling.unwrap() - 4.52).abs() < 0.01);

        // Set experimental chemical shift, partial RMSD should be within 0.1 of 43.0
        residue.set_exp_shift(135.499);
        assert!((residue.rmsd.unwrap() - 43.0).abs() < 0.1);

        // Set experimental dipolar coupling, weighted RMSD should be within 0.1 of 37.8
        residue.set_exp_dipolar_coupling(1.348);
        assert!((residue.rmsd.unwrap() - 37.8).abs() < 0.1);

        // Unset experimental chemical shift, partial RMSD should be within 0.1 of 31.7
        residue.unset_exp_shift();
        assert!((residue.rmsd.unwrap() - 31.7).abs() < 0.1);

        // Unset experimental dipolar coupling, RMSD should be empty again
        residue.unset_exp_dipolar_coupling();
        assert!(residue.rmsd.is_none());
    }
}
