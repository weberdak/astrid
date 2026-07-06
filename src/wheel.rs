use crate::matrix_a::calc_matrix_a;
use crate::params::Params;
use crate::residue::{self, Residue, ResidueData};
use nalgebra::Matrix3;

/// A wheel of residues with experimental and simulated peaks
/// 
/// Methods provide for dynamic and efficient recalculation of simulated peaks and RMSD 
/// when parameters are changed one at a time.
pub struct Wheel {
    /// Parameters for the simulation
    params: Params,
    /// Matrix A for the wheel
    matrix_a: Matrix3<f64>,
    /// Number of residues in the wheel
    num_residues: usize,
    /// Residue number of the reference residue in the wheel
    residue_number_ref: usize,
    /// Residue number of the first residue in the wheel
    residue_number_first: usize,
    /// Residues in the wheel
    pub residues: Vec<Residue>,
    /// Mean RMSD over all calculated and simulated peaks
    pub rmsd: f64,
}

impl Wheel {
    /// Create a new Wheel with the given parameters and number of residues
    pub fn new(
        params: Params,
        num_residues: usize,
        residue_number_ref: usize,
        residue_number_first: usize,
    ) -> Self {
        let matrix_a = calc_matrix_a(&params);
        let offset: isize = residue_number_ref as isize - residue_number_first as isize;
        let residues: Vec<Residue> = (0..num_residues)
            .map(|i| Residue::new(i as isize - offset, &params))
            .collect();
        let mut wheel = Self {
            params,
            matrix_a,
            num_residues,
            residue_number_ref,
            residue_number_first,
            residues,
            rmsd: 0.0,
        };
        wheel.update_sim_peaks();
        wheel
    }

    /// Create a wheel pre-populated with experimental data
    pub fn from_data(
        params: Params,
        residue_number_ref: usize,
        residue_number_first: usize,
        residue_data: Vec<ResidueData>,
    ) -> Self {
        let matrix_a = calc_matrix_a(&params);
        let residues: Vec<Residue> = residue_data
            .into_iter()
            .map(|data| Residue::from_data(data, &params, &matrix_a))
            .collect();
        let num_residues = residues.len();
        let mut wheel = Self {
            params,
            matrix_a,
            num_residues,
            residue_number_ref,
            residue_number_first,
            residues,
            rmsd: 0.0,
        };
        wheel.update_rmsd();
        wheel
    }

    /// Set or unset label of specific residue in the wheel
    pub fn update_label(&mut self, residue_number: usize, label: Option<[u8; 1]>) {
        self.validate_residue_number(residue_number);
        let residue_index = residue_number - self.residue_number_first;
        self.residues[residue_index].label = label;
    }

    /// Set or unset experimental shift of specific residue in the wheel
    pub fn update_exp_shift(&mut self, residue_number: usize, exp_shift: Option<f64>) {
        self.validate_residue_number(residue_number);
        let residue_index = residue_number - self.residue_number_first;
        self.residues[residue_index].update_exp_shift(exp_shift);
        self.update_rmsd();
    }

    /// Set or unset experimental dipolar coupling of specific residue in the wheel
    pub fn update_exp_dipolar_coupling(
        &mut self,
        residue_number: usize,
        exp_dipolar_coupling: Option<f64>,
    ) {
        self.validate_residue_number(residue_number);
        let residue_index = residue_number - self.residue_number_first;
        self.residues[residue_index].update_exp_dipolar_coupling(exp_dipolar_coupling);
        self.update_rmsd();
    }

    /// Update rotation of reference residue and recalculate simulated peaks
    pub fn update_rotation(&mut self, rotation_deg: f64) {
        self.params.rotation = rotation_deg.to_radians();
        for residue in &mut self.residues {
            residue.update_rotation(&self.params, &self.matrix_a);
        }
        self.update_rmsd();
    }

    /// Update tilt angle of helix and recalculate simulated peaks
    pub fn update_tilt(&mut self, tilt_deg: f64) {
        // Fold angle to 0-90 degrees range
        let normalised = tilt_deg.rem_euclid(360.0);
        let mirrored = if normalised > 180.0 {
            360.0 - normalised
        } else {
            normalised
        };
        let folded = if mirrored > 90.0 {
            180.0 - mirrored
        } else {
            mirrored
        };
        self.params.tilt = folded.to_radians();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update order parameter of helix and recalculate simulated peaks
    pub fn update_order_parameter(&mut self, order_parameter: f64) {
        if !(0.0..=90.0).contains(&order_parameter) {
            panic!(
                "Order parameter {} is out of range (must be between 0 and 90)",
                order_parameter
            );
        }
        self.params.update_order_parameter(order_parameter);
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update flip angle of helix and recalculate simulated peaks
    pub fn update_flip_angle(&mut self, flip_angle_deg: f64) {
        if !(0.0..=90.0).contains(&flip_angle_deg) {
            panic!(
                "Flip angle {} is out of range (must be between 0 and 90)",
                flip_angle_deg
            );
        }
        self.params.update_flip_angle(flip_angle_deg);
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update reference residue number
    pub fn update_residue_number_ref(&mut self, residue_number_ref: usize) {
        self.validate_residue_number(residue_number_ref);
        self.residue_number_ref = residue_number_ref;
        let offset: isize =
            residue_number_ref as isize - self.residue_number_first as isize;
        for (i, residue) in self.residues.iter_mut().enumerate() {
            residue.index = i as isize - offset;
        }
        self.params.rotation = self.residues[offset as usize].rotation;
    }

    /// Update phi dihedral angle of helix and recalculate simulated peaks
    pub fn update_phi(&mut self, phi_deg: f64) {
        if !(-180.0..=180.0).contains(&phi_deg) {
            panic!(
                "Phi angle {} is out of range (must be between -180 and 180)",
                phi_deg
            );
        }
        self.params.phi = phi_deg.to_radians();
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update psi dihedral angle of helix and recalculate simulated peaks
    pub fn update_psi(&mut self, psi_deg: f64) {
        if !(-180.0..=180.0).contains(&psi_deg) {
            panic!(
                "Psi angle {} is out of range (must be between -180 and 180)",
                psi_deg
            );
        }
        self.params.psi = psi_deg.to_radians();
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update beta rotation angle and recalculate simulated peaks
    pub fn update_beta(&mut self, beta_deg: f64) {
        self.params.beta = beta_deg.to_radians();
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update CA-C-N bond angle and recalculate simulated peaks
    pub fn update_ca_c_n(&mut self, ca_c_n_deg: f64) {
        if !(0.0..=180.0).contains(&ca_c_n_deg) {
            panic!(
                "CA-C-N bond angle {} is out of range (must be between 0 and 180)",
                ca_c_n_deg
            );
        }
        self.params.ca_c_n = ca_c_n_deg.to_radians();
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update C-N-CA bond angle and recalculate simulated peaks
    pub fn update_c_n_ca(&mut self, c_n_ca_deg: f64) {
        if !(0.0..=180.0).contains(&c_n_ca_deg) {
            panic!(
                "C-N-CA bond angle {} is out of range (must be between 0 and 180)",
                c_n_ca_deg
            );
        }
        self.params.c_n_ca = c_n_ca_deg.to_radians();
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update N-CA-C bond angle and recalculate simulated peaks
    pub fn update_n_ca_c(&mut self, n_ca_c_deg: f64) {
        if !(0.0..=180.0).contains(&n_ca_c_deg) {
            panic!(
                "N-CA-C bond angle {} is out of range (must be between 0 and 180)",
                n_ca_c_deg
            );
        }
        self.params.n_ca_c = n_ca_c_deg.to_radians();
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update CA-N-H bond angle and recalculate simulated peaks
    pub fn update_ca_n_h(&mut self, ca_n_h_deg: f64) {
        if !(0.0..=180.0).contains(&ca_n_h_deg) {
            panic!(
                "CA-N-H bond angle {} is out of range (must be between 0 and 180)",
                ca_n_h_deg
            );
        }
        self.params.ca_n_h = ca_n_h_deg.to_radians();
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update the CA-C bond length and recalculate simulated peaks
    pub fn update_ca_c(&mut self, ca_c_length: f64) {
        if ca_c_length < 0.0 {
            panic!(
                "CA-C bond length {} is negative (must be non-negative)",
                ca_c_length
            );
        }
        self.params.ca_c = ca_c_length;
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update the C-N bond length and recalculate simulated peaks
    pub fn update_c_n(&mut self, c_n_length: f64) {
        if c_n_length < 0.0 {
            panic!(
                "C-N bond length {} is negative (must be non-negative)",
                c_n_length
            );
        }
        self.params.c_n = c_n_length;
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update the N-CA bond length and recalculate simulated peaks
    pub fn update_n_ca(&mut self, n_ca_length: f64) {
        if n_ca_length < 0.0 {
            panic!(
                "N-CA bond length {} is negative (must be non-negative)",
                n_ca_length
            );
        }
        self.params.n_ca = n_ca_length;
        self.update_matrix_a();
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update the δxx shift tensor component and recalculate simulated peaks
    pub fn update_delta_xx(&mut self, delta_xx: f64) {
        self.params.update_delta_xx(delta_xx);
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update the δyy shift tensor component and recalculate simulated peaks
    pub fn update_delta_yy(&mut self, delta_yy: f64) {
        self.params.update_delta_yy(delta_yy);
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update the δzz shift tensor component and recalculate simulated peaks
    pub fn update_delta_zz(&mut self, delta_zz: f64) {
        self.params.update_delta_zz(delta_zz);
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Update dipolar coupling constant and recalculate simulated peaks
    pub fn update_coupling_constant(&mut self, coupling: f64) {
        if (coupling < 0.0) {
            panic!(
                "Dipolar coupling constant {} is negative (must be non-negative)",
                coupling
            );
        }
        self.params.coupling = coupling;
        self.update_sim_peaks();
        self.update_rmsd();
    }

    /// Print label+number, rotation, experimental and simulated peaks for all residues in the wheel
    pub fn as_csv(&self, sig_figs: usize) -> String {
        let mut csv_data = String::from(
            "Residue,Rotation,Exp Shift,Sim Shift,Exp Coupling,Sim Coupling,RMSD\n",
        );
        for residue in &self.residues {
            let format_sig_figs = |opt: Option<f64>| -> String {
                match opt {
                    Some(val) => {
                        if val == 0.0 {
                            "0".to_string()
                        } else {
                            let magnitude = val.abs().log10().floor() as i32;
                            let decimals =
                                (sig_figs as i32 - magnitude - 1).max(0) as usize;
                            format!("{:.prec$}", val, prec = decimals)
                        }
                    }
                    None => "None".to_string(),
                }
            };
            let residue_number = residue.index + self.residue_number_ref as isize;
            let residue_label = match residue.label {
                Some(label) => format!("{}{}", label[0] as char, residue_number),
                None => format!("{}", residue_number),
            };
            csv_data.push_str(&format!(
                "{},{},{},{},{},{},{}\n",
                residue_label,
                format_sig_figs(Some(residue.rotation.to_degrees())),
                format_sig_figs(residue.exp_shift),
                format_sig_figs(residue.sim_shift),
                format_sig_figs(residue.exp_dipolar_coupling),
                format_sig_figs(residue.sim_dipolar_coupling),
                format_sig_figs(residue.rmsd),
            ));
        }
        csv_data
    }

    /// Update Matrix A (peptide bond and angle params change)
    fn update_matrix_a(&mut self) {
        self.matrix_a = calc_matrix_a(&self.params);
    }

    /// Validate that the residue number is within bounds for the wheel
    fn validate_residue_number(&self, residue_number: usize) {
        if residue_number < self.residue_number_first {
            panic!(
                "Residue number {} is too low (minimum is {})",
                residue_number, self.residue_number_first
            );
        }

        let max_residue = self.residue_number_first + self.num_residues - 1;
        if residue_number > max_residue {
            panic!(
                "Residue number {} is too high (maximum is {})",
                residue_number, max_residue
            );
        }
    }

    /// Update the simulated peaks for all residues in the wheel
    fn update_sim_peaks(&mut self) {
        for residue in &mut self.residues {
            residue.calc_peak(&self.params, &self.matrix_a);
        }
    }

    /// Update the RMSD for the wheel
    fn update_rmsd(&mut self) {
        let (total, count) = self
            .residues
            .iter()
            .filter(|residue| !residue.ignore_rmsd)
            .filter_map(|residue| residue.rmsd)
            .fold((0.0, 0usize), |(total, count), rmsd| {
                (total + rmsd, count + 1)
            });

        self.rmsd = if count > 0 { total / count as f64 } else { 0.0 };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;

    #[test]
    /// Update label of specific residue in the wheel
    fn test_update_residue_label() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_label(2, Some([b'A']));
        assert_eq!(wheel.residues[1].label, Some([b'A']));
        wheel.update_label(2, None);
        assert_eq!(wheel.residues[1].label, None);
    }

    #[test]
    /// Raise panic when residue number is too low
    fn test_update_label_low_residue_number() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_label(0, Some([b'A']));
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Raise panic when residue number is too high
    fn test_update_label_high_residue_number() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_label(6, Some([b'A']));
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Update experimental shift of specific residue in the wheel
    fn test_update_exp_shift() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(2, Some(200.0));
        assert_eq!(wheel.residues[1].exp_shift, Some(200.0));
        assert_eq!(wheel.rmsd, wheel.residues[1].rmsd.unwrap());
    }

    #[test]
    /// Raise panic when residue number is too low
    fn test_update_exp_shift_low_residue_number() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_exp_shift(0, Some(200.0));
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Raise panic when residue number is too high
    fn test_update_exp_shift_high_residue_number() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_exp_shift(6, Some(200.0));
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Update experimental dipolar coupling of specific residue in the wheel
    fn test_update_exp_dipolar_coupling() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_dipolar_coupling(2, Some(4.0));
        assert_eq!(wheel.residues[1].exp_dipolar_coupling, Some(4.0));
        assert_eq!(wheel.rmsd, wheel.residues[1].rmsd.unwrap());
    }

    #[test]
    /// Raise panic when residue number is too low
    fn test_update_exp_dipolar_coupling_low_residue_number() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_exp_dipolar_coupling(0, Some(4.0));
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Raise panic when residue number is too high
    fn test_update_exp_dipolar_coupling_high_residue_number() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_exp_dipolar_coupling(6, Some(4.0));
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Update rotation of reference residue and recalculate simulated peaks
    fn test_update_rotation_expected_direction() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let residue_1_rotation = wheel.residues[0].rotation;
        let residue_1_sim_shift = wheel.residues[0].sim_shift;
        let residue_1_sim_coupling = wheel.residues[0].sim_dipolar_coupling;
        wheel.update_rotation(100.0); // degrees to next residue on the helix
        // Second residue should match first residue prior to rotation
        assert_eq!(wheel.residues[1].rotation, residue_1_rotation);
        assert_eq!(wheel.residues[1].sim_shift, residue_1_sim_shift);
        assert_eq!(
            wheel.residues[1].sim_dipolar_coupling,
            residue_1_sim_coupling
        );
    }

    #[test]
    /// Update rotation of reference residue and recalculate simulated peaks
    fn test_update_rotation_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.rotation.to_degrees() - 0.0).abs() < 0.1);
        wheel.update_rotation(100.0);
        assert!((wheel.params.rotation.to_degrees() - 100.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 224.23).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 7.80).abs() < 0.1);
        assert!((wheel.rmsd - 20.65).abs() < 0.1);
    }

    #[test]
    /// Residues re-indexed without changing simulated peak values
    fn test_update_residue_number_ref() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let residue_3_rotation = wheel.residues[2].rotation;
        let residue_1_shift_before = wheel.residues[0].sim_shift;
        let residue_1_coupling_before = wheel.residues[0].sim_dipolar_coupling;
        assert_eq!(wheel.residues[0].index, 0);
        wheel.update_residue_number_ref(3);
        assert_eq!(wheel.residues[0].index, -2);
        assert_eq!(wheel.residues[0].sim_shift, residue_1_shift_before);
        assert_eq!(
            wheel.residues[0].sim_dipolar_coupling,
            residue_1_coupling_before
        );
        assert_eq!(wheel.params.rotation, residue_3_rotation);
    }

    #[test]
    /// Tilt angle of helix updated, peaks and RMSD properly recalculated
    fn test_update_tilt_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.tilt.to_degrees() - 20.0).abs() < 0.1);
        wheel.update_tilt(30.0);
        assert!((wheel.params.tilt.to_degrees() - 30.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 182.96).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 9.64).abs() < 0.1);
        assert!((wheel.rmsd - 15.79).abs() < 0.1);
    }

    #[test]
    /// Tilt angle of helix updated with correct folding
    fn test_update_tilt_folding() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_tilt(100.0);
        assert!((wheel.params.tilt.to_degrees() - 80.0).abs() < 0.1);
        wheel.update_tilt(200.0);
        assert!((wheel.params.tilt.to_degrees() - 20.0).abs() < 0.1);
        wheel.update_tilt(280.0);
        assert!((wheel.params.tilt.to_degrees() - 80.0).abs() < 0.1);
        wheel.update_tilt(350.0);
        assert!((wheel.params.tilt.to_degrees() - 10.0).abs() < 0.1);
    }

    #[test]
    /// Order parameter of helix updated, peaks and RMSD properly recalculated
    fn test_update_order_parameter_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.order - 1.0).abs() < 0.1);
        wheel.update_order_parameter(0.8);
        assert!((wheel.params.order - 0.8).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 187.239).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 8.47).abs() < 0.1);
        assert!((wheel.rmsd - 16.56).abs() < 0.1);
    }

    #[test]
    /// Order parameter of helix updated with out-of-range values, panics as expected
    fn test_update_order_parameter_out_of_range() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_order_parameter(-0.1);
        }));
        assert!(result.is_err());
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_order_parameter(90.1);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Flip angle of helix updated, peaks and RMSD properly recalculated
    fn test_update_flip_angle_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.flip.to_degrees() - 0.0).abs() < 0.1);
        wheel.update_flip_angle(90.0);
        assert!((wheel.params.flip.to_degrees() - 90.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 81.55).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 5.30).abs() < 0.1);
        assert!((wheel.rmsd - 93.4).abs() < 0.1);
    }

    #[test]
    /// Flip angle of helix updated with out-of-range values, panics as expected
    fn test_update_flip_angle_out_of_range() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_flip_angle(-0.1);
        }));
        assert!(result.is_err());
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_flip_angle(90.1);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Phi angle of helix updated, peaks and RMSD properly recalculated
    fn test_update_phi_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.phi.to_degrees() - (-63.0)).abs() < 0.1);
        wheel.update_phi(-70.0);
        assert!((wheel.params.phi.to_degrees() - (-70.0)).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 198.69).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.44).abs() < 0.1);
        assert!((wheel.rmsd - 5.45).abs() < 0.1);
    }

    #[test]
    /// Phi angle of helix updated with out-of-range values, panics as expected
    fn test_update_phi_out_of_range() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_phi(-180.1);
        }));
        assert!(result.is_err());
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_phi(180.1);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Psi angle of helix updated, peaks and RMSD properly recalculated
    fn test_update_psi_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.psi.to_degrees() - (-42.0)).abs() < 0.1);
        wheel.update_psi(-50.0);
        assert!((wheel.params.psi.to_degrees() - (-50.0)).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 196.95).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.34).abs() < 0.1);
        assert!((wheel.rmsd - 6.18).abs() < 0.1);
    }

    #[test]
    /// Psi angle of helix updated with out-of-range values, panics as expected
    fn test_update_psi_out_of_range() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_psi(-180.1);
        }));
        assert!(result.is_err());
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_psi(180.1);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// Beta angle updated, peaks and RMSD properly recalculated
    fn test_update_beta_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.beta.to_degrees() - 17.0).abs() < 0.1);
        wheel.update_beta(25.0);
        assert!((wheel.params.beta.to_degrees() - 25.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 184.63).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.59).abs() < 0.1);
        assert!((wheel.rmsd - 15.00).abs() < 0.1);
    }

    #[test]
    /// CA-C-N bond angle updated, peaks and RMSD properly recalculated
    fn test_update_ca_c_n_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.ca_c_n.to_degrees() - 117.5).abs() < 0.1);
        wheel.update_ca_c_n(125.0);
        assert!((wheel.params.ca_c_n.to_degrees() - 125.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 210.41).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.62).abs() < 0.1);
        assert!((wheel.rmsd - 5.87).abs() < 0.1);
    }

    #[test]
    /// CA-C-N bond angle updated with out-of-range values, panics as expected
    fn test_update_ca_c_n_out_of_range() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_ca_c_n(-0.1);
        }));
        assert!(result.is_err());
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_ca_c_n(180.1);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// C-N-CA bond angle updated, peaks and RMSD properly recalculated
    fn test_update_c_n_ca_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.c_n_ca.to_degrees() - 124.0).abs() < 0.1);
        wheel.update_c_n_ca(145.0);
        assert!((wheel.params.c_n_ca.to_degrees() - 145.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 173.75).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 9.14).abs() < 0.1);
        assert!((wheel.rmsd - 22.91).abs() < 0.1);
    }

    #[test]
    /// C-N-CA bond angle updated with out-of-range values, panics as expected
    fn test_update_c_n_ca_out_of_range() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_c_n_ca(-0.1);
        }));
        assert!(result.is_err());
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_c_n_ca(180.1);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// N-CA-C bond angle updated, peaks and RMSD properly recalculated
    fn test_update_n_ca_c_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.n_ca_c.to_degrees() - 107.4).abs() < 0.1);
        wheel.update_n_ca_c(115.0);
        assert!((wheel.params.n_ca_c.to_degrees() - 115.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 211.06).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.72).abs() < 0.1);
        assert!((wheel.rmsd - 6.64).abs() < 0.1);
    }

    #[test]
    /// N-CA-C bond angle updated with out-of-range values, panics as expected
    fn test_update_n_ca_c_out_of_range() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_n_ca_c(-0.1);
        }));
        assert!(result.is_err());
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_n_ca_c(180.1);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// CA-N-H bond angle updated, peaks and RMSD properly recalculated
    fn test_update_ca_n_h_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.ca_n_h.to_degrees() - 116.0).abs() < 0.1);
        wheel.update_ca_n_h(125.0);
        assert!((wheel.params.ca_n_h.to_degrees() - 125.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 182.01).abs() < 0.1);
        assert!((wheel.residues[0].sim_dipolar_coupling.unwrap() - 9.75).abs() < 0.1);
        assert!((wheel.rmsd - 16.35).abs() < 0.1);
    }

    #[test]
    /// CA-N-H bond angle updated with out-of-range values, panics as expected
    fn test_update_ca_n_h_out_of_range() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_ca_n_h(-0.1);
        }));
        assert!(result.is_err());
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_ca_n_h(180.1);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// CA-C bond length updated, peaks and RMSD properly recalculated
    fn test_update_ca_c_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.ca_c - 1.52).abs() < 0.1);
        wheel.update_ca_c(1.7);
        assert!((wheel.params.ca_c - 1.7).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 203.418).abs() < 0.001);
        assert!(
            (wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.5902).abs() < 0.0001
        );
        assert!((wheel.rmsd - 4.3209).abs() < 0.0001);
    }

    #[test]
    /// CA-C bond length updated with negative value, panics as expected
    fn test_update_ca_c_negative_value() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_ca_c(-1.0);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// C-N bond length updated, peaks and RMSD properly recalculated
    fn test_update_c_n_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.c_n - 1.35).abs() < 0.1);
        wheel.update_c_n(1.9);
        assert!((wheel.params.c_n - 1.9).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 203.919).abs() < 0.001);
        assert!(
            (wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.5989).abs() < 0.0001
        );
        assert!((wheel.rmsd - 4.3034).abs() < 0.0001);
    }

    #[test]
    /// C-N bond length updated with negative value, panics as expected
    fn test_update_c_n_negative_value() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_c_n(-1.0);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// N-CA bond length updated, peaks and RMSD properly recalculated
    fn test_update_n_ca_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.n_ca - 1.45).abs() < 0.1);
        wheel.update_n_ca(1.9);
        assert!((wheel.params.n_ca - 1.9).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 203.367).abs() < 0.001);
        assert!(
            (wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.5893).abs() < 0.0001
        );
        assert!((wheel.rmsd - 4.3239).abs() < 0.0001);
    }

    #[test]
    /// N-CA bond length updated with negative value, panics as expected
    fn test_update_n_ca_negative_value() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_n_ca(-1.0);
        }));
        assert!(result.is_err());
    }

    #[test]
    /// δxx shift tensor component updated, peaks and RMSD properly recalculated
    fn test_update_delta_xx_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        let shift_iso_before = wheel.params.shift_iso;
        assert!((wheel.params.delta_xx - 57.3).abs() < 0.1);
        wheel.update_delta_xx(20.0);
        assert!((wheel.params.shift_iso - shift_iso_before).abs() > 0.1);
        assert!((wheel.params.delta_xx - 20.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 198.14).abs() < 0.1);
        assert!(
            (wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.5917).abs() < 0.0001
        );
        assert!((wheel.rmsd - 6.4057).abs() < 0.1);
    }

    #[test]
    /// δyy shift tensor component updated, peaks and RMSD properly recalculated
    fn test_update_delta_yy_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        let shift_iso_before = wheel.params.shift_iso;
        assert!((wheel.params.delta_yy - 81.2).abs() < 0.1);
        wheel.update_delta_yy(100.0);
        assert!((wheel.params.shift_iso - shift_iso_before).abs() > 0.1);
        assert!((wheel.params.delta_yy - 100.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 203.507).abs() < 0.001);
        assert!(
            (wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.5917).abs() < 0.0001
        );
        assert!((wheel.rmsd - 4.3152).abs() < 0.0001);
    }

    #[test]
    /// δzz shift tensor component updated, peaks and RMSD properly recalculated
    fn test_update_delta_zz_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        let shift_iso_before = wheel.params.shift_iso;
        assert!((wheel.params.delta_zz - 228.1).abs() < 0.1);
        wheel.update_delta_zz(250.0);
        assert!((wheel.params.shift_iso - shift_iso_before).abs() > 0.1);
        assert!((wheel.params.delta_zz - 250.0).abs() < 0.1);
        assert!((wheel.residues[0].sim_shift.unwrap() - 222.243).abs() < 0.1);
        assert!(
            (wheel.residues[0].sim_dipolar_coupling.unwrap() - 10.5917).abs() < 0.0001
        );
        assert!((wheel.rmsd - 12.89).abs() < 0.1);
    }

    #[test]
    /// Dipolar coupling constant updated, peaks and RMSD properly recalculated
    fn test_update_dipolar_coupling_constant_expected_values() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        wheel.update_exp_shift(1, Some(205.0));
        wheel.update_exp_dipolar_coupling(1, Some(10.0));
        assert!((wheel.params.coupling - 10.735).abs() < 0.001);
        wheel.update_coupling_constant(15.0);
        assert!((wheel.params.coupling - 15.0).abs() < 0.0001);
        assert!((wheel.residues[0].sim_shift.unwrap() - 203.499).abs() < 0.001);
        assert!(
            (wheel.residues[0].sim_dipolar_coupling.unwrap() - 14.7998).abs() < 0.0001
        );
        assert!((wheel.rmsd - 33.9563).abs() < 0.0001);
    }

    #[test]
    /// Dipolar coupling constant updated with negative value, panics as expected
    fn test_update_dipolar_coupling_constant_negative_value() {
        let params = Params::new();
        let mut wheel = Wheel::new(params, 5, 1, 1);
        let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
            wheel.update_coupling_constant(-1.0);
        }));
        assert!(result.is_err());
    }
}
