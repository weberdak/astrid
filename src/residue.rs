use crate::paf::calc_paf_coords;
use crate::params::Params;
use crate::peak::{calc_dipolar_coupling, calc_shift};
use nalgebra::Matrix3;
use std::f64::consts::PI;

/// Residue initialisation data from user input
pub struct ResidueData {
    pub index: isize,
    pub label: Option<[u8; 1]>,
    pub exp_shift: Option<f64>,
    pub exp_dipolar_coupling: Option<f64>,
    pub ignore_rmsd: bool,
}

impl ResidueData {
    /// Create a new ResidueData with user-defined values.
    pub fn new(
        index: isize,
        label: Option<[u8; 1]>,
        exp_shift: Option<f64>,
        exp_dipolar_coupling: Option<f64>,
        ignore_rmsd: bool,
    ) -> Self {
        Self {
            index,
            label,
            exp_shift,
            exp_dipolar_coupling,
            ignore_rmsd,
        }
    }
}

/// Stores labelling and experimental data for a residue in a protein sequence, and
/// calculates simulated peaks and RMSD.
pub struct Residue {
    /// Index relative to reference residue in the sequence
    pub index: isize,
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
    /// Ignore RMSD for this residue in the overall wheel RMSD calculation
    pub ignore_rmsd: bool,
}

impl Residue {
    /// Create a new Residue instance with default values
    pub fn new(index: isize, params: &Params) -> Self {
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
            ignore_rmsd: false,
        }
    }

    /// Create a new Residue instance from ResidueData and calculate simulated peaks
    pub fn from_data(
        data: ResidueData,
        params: &Params,
        matrix_a: &Matrix3<f64>,
    ) -> Self {
        let mut residue = Self::new(data.index, params);
        residue.label = data.label;
        residue.exp_shift = data.exp_shift;
        residue.exp_dipolar_coupling = data.exp_dipolar_coupling;
        residue.ignore_rmsd = data.ignore_rmsd;
        residue.calc_peak(params, matrix_a);
        residue
    }

    /// Update the simulated shift and dipolar coupling for the residue
    pub fn calc_peak(&mut self, params: &Params, matrix_a: &Matrix3<f64>) {
        let paf_coords = calc_paf_coords(self.rotation, params, matrix_a);
        self.sim_shift = Some(calc_shift(params, &paf_coords));
        self.sim_dipolar_coupling =
            Some(calc_dipolar_coupling(params, &paf_coords).abs());
        self.calc_rmsd();
    }

    /// Update the rotation of the reference residue and recalculate simulated peaks
    pub fn update_rotation(&mut self, params: &Params, matrix_a: &Matrix3<f64>) {
        self.rotation = (params.rotation - self.index as f64 * 100.0_f64.to_radians())
            .rem_euclid(2.0 * PI);
        self.calc_peak(params, matrix_a);
    }

    /// Update the experimental shift for the residue and recalculate RMSD
    pub fn update_exp_shift(&mut self, exp_shift: Option<f64>) {
        self.exp_shift = exp_shift;
        self.calc_rmsd();
    }

    /// Update the experimental dipolar coupling for the residue and recalculate RMSD
    pub fn update_exp_dipolar_coupling(&mut self, exp_dipolar_coupling: Option<f64>) {
        self.exp_dipolar_coupling = exp_dipolar_coupling;
        self.calc_rmsd();
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
    ///
    /// This RMSD calculation differs from the original Python code
    /// (https://github.com/weberdak/pisa.py/) and has been corrected so that the
    /// dipolar coupling deviation is scaled rather than the root of squared dipolar
    /// coupling deviation.
    fn calc_rmsd(&mut self) {
        // Factor used to match magnitudes of dipolar coupling and chemical shift
        // deviations.
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

    /// I3 of Sarcolipin (SLN) from Weber et al. (2020)
    /// https://github.com/weberdak/pisa.py/tree/master/examples/sarcolipin/sln_explore_log.dat
    fn generate_test_params() -> Params {
        let mut params = Params::new();
        params.update_flip_angle(90.0);
        params.update_order_parameter(0.9);
        params.tilt = 24.6_f64.to_radians();
        params.rotation = 46.0_f64.to_radians();
        params
    }

    fn generate_test_matrix_a() -> Matrix3<f64> {
        let params = generate_test_params();
        calc_matrix_a(&params)
    }

    fn generate_test_residue() -> Residue {
        let params = generate_test_params();
        Residue::new(-3, &params)
    }

    #[test]
    /// Rotation should always be between 0 and 2π radians
    fn test_rotation_within_bounds_at_init() {
        let mut params = generate_test_params();
        params.rotation = 1000.0_f64.to_radians();
        let residue = Residue::new(-3, &params);
        assert!(residue.rotation >= 0.0 && residue.rotation < 2.0 * PI);
        params.rotation = -1000.0_f64.to_radians();
        let residue = Residue::new(-3, &params);
        assert!(residue.rotation >= 0.0 && residue.rotation < 2.0 * PI);
    }

    #[test]
    /// Rotation update applied correctly
    fn test_rotation_update_is_expected_value() {
        let matrix_a = generate_test_matrix_a();
        let mut params = generate_test_params();
        params.rotation = 0.0_f64.to_radians();
        // Second residue is 100 deg less than the reference residue
        let mut residue = Residue::new(1, &params);
        assert!((residue.rotation - 260.0_f64.to_radians()).abs() < 1e-6);
        params.rotation = 100.0_f64.to_radians();
        residue.update_rotation(&params, &matrix_a); // Should reduce to 0 deg
        assert!((residue.rotation - 0.0).abs() < 1e-6);
    }

    #[test]
    /// Rotation should always be between 0 and 2π radians
    fn test_rotation_within_bounds_after_update() {
        let matrix_a = generate_test_matrix_a();
        let mut params = generate_test_params();
        params.rotation = 1000.0_f64.to_radians();
        let mut residue = Residue::new(-3, &params);
        residue.update_rotation(&params, &matrix_a);
        assert!(residue.rotation >= 0.0 && residue.rotation < 2.0 * PI);
        params.rotation = -1000.0_f64.to_radians();
        let mut residue = Residue::new(-3, &params);
        residue.update_rotation(&params, &matrix_a);
        assert!(residue.rotation >= 0.0 && residue.rotation < 2.0 * PI);
    }

    #[test]
    /// RMSD should stay None if no experimental values are set
    fn test_rmsd_none_after_calc_peak_without_exp_values() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        assert!(residue.rmsd.is_none());
    }

    #[test]
    /// Simulated shift matches expected value
    fn test_sim_shift_is_expected_value() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        assert!((residue.sim_shift.unwrap() - 92.52).abs() < 0.1);
    }

    #[test]
    /// Simulated dipolar coupling matches expected value and converted to positive
    fn test_sim_dipolar_coupling_is_expected_value() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        assert!((residue.sim_dipolar_coupling.unwrap() - 4.52).abs() < 0.01);
    }

    #[test]
    /// Partial RMSD is calculated when only experimental shift is set
    fn test_partial_rmsd_with_only_exp_shift() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        residue.update_exp_shift(Some(135.499));
        assert!((residue.rmsd.unwrap() - 43.0).abs() < 0.1);
    }

    #[test]
    /// Partial RMSD is calculated correctly when only dipolar coupling is set
    fn test_partial_rmsd_with_only_exp_dipolar_coupling() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        residue.update_exp_dipolar_coupling(Some(1.348));
        assert!((residue.rmsd.unwrap() - 31.7).abs() < 0.1);
    }

    #[test]
    /// Full RMSD is calculated correctly when both experimental shift and dipolar coupling are set
    fn test_rmsd_with_shift_and_dipolar_coupling() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        residue.update_exp_shift(Some(135.499));
        residue.update_exp_dipolar_coupling(Some(1.348));
        assert!((residue.rmsd.unwrap() - 37.8).abs() < 0.1);
    }

    #[test]
    /// RMSD returns to partial when shift is unset
    fn test_rmsd_returns_to_partial_when_shift_unset() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        residue.update_exp_shift(Some(135.499));
        residue.update_exp_dipolar_coupling(Some(1.348));
        residue.update_exp_shift(None);
        assert!((residue.rmsd.unwrap() - 31.7).abs() < 0.1);
    }

    #[test]
    /// RMSD returns to partial when dipolar coupling is unset
    fn test_rmsd_returns_to_partial_when_dipolar_coupling_unset() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        residue.update_exp_shift(Some(135.499));
        residue.update_exp_dipolar_coupling(Some(1.348));
        residue.update_exp_dipolar_coupling(None);
        assert!((residue.rmsd.unwrap() - 43.0).abs() < 0.1);
    }

    #[test]
    /// RMSD returns to None when both experimental values are unset
    fn test_rmsd_returns_to_none_when_both_unset() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        residue.update_exp_shift(Some(135.499));
        residue.update_exp_dipolar_coupling(Some(1.348));
        assert!((residue.rmsd.unwrap() - 37.8).abs() < 0.1);
        residue.update_exp_shift(None);
        residue.update_exp_dipolar_coupling(None);
        assert!(residue.rmsd.is_none());
    }

    #[test]
    fn rmsd_calculated_despite_ignore_rmsd() {
        let mut residue = generate_test_residue();
        residue.calc_peak(&generate_test_params(), &generate_test_matrix_a());
        residue.ignore_rmsd = true;
        residue.update_exp_shift(Some(135.499));
        residue.update_exp_dipolar_coupling(Some(1.348));
        assert!((residue.rmsd.unwrap() - 37.8).abs() < 0.1);
    }
}
