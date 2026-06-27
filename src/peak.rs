use crate::params::Params;
use nalgebra::Vector3;

/// Calculate the oriented 15N chemical shift in ppm
pub fn calc_shift(params: &Params, paf_coords: &Vector3<f64>) -> f64 {
    // Value without factoring in order or flip angle.
    let raw = params.delta_xx * paf_coords.x.powi(2)
        + params.delta_yy * paf_coords.y.powi(2)
        + params.delta_zz * paf_coords.z.powi(2);
    (raw - params.shift_iso) * params.scalar + params.shift_iso
}

/// Calculate the oriented 15N-1H dipolar coupling in kHz
pub fn calc_dipolar_coupling(params: &Params, paf_coords: &Vector3<f64>) -> f64 {
    params.coupling
        * 0.5
        * (3.0
            * (params.beta.sin() * paf_coords.x + params.beta.cos() * paf_coords.z)
                .powi(2)
            - 1.0)
        * params.scalar
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix_a::calc_matrix_a;
    use crate::paf::calc_paf_coords;
    use crate::params::Params;

    #[test]
    fn test_calc_paf_coords() {
        let params = Params::new();
        let matrix_a = calc_matrix_a(&params);

        // [rotation, shift, coupling] for default parameters
        let expected = [
            [0.0, 203.5, 10.59],
            [4.53778, 184.8, 7.06],
            [2.79248, 208.7, 5.67],
            [1.04718, 223.5, 9.55],
        ];

        let tol = 1e-1;

        for row in expected {
            let (rotation, expected_shift, expected_coupling) =
                (row[0], row[1], row[2]);

            let paf_coords = calc_paf_coords(rotation, &params, &matrix_a);
            let shift = calc_shift(&params, &paf_coords);
            let coupling = calc_dipolar_coupling(&params, &paf_coords);

            assert!(
                (shift - expected_shift).abs() < tol,
                "shift mismatch at rotation {}: got {}, expected {}",
                rotation,
                shift,
                expected_shift
            );

            assert!(
                (coupling - expected_coupling).abs() < tol,
                "coupling mismatch at rotation {}: got {}, expected {}",
                rotation,
                coupling,
                expected_coupling
            );
        }
    }

    #[test]
    fn test_calc_paf_coords_flipped() {
        let mut params = Params::new();
        params.update_flip_angle(90.0);
        params.update_order_parameter(0.8);
        let matrix_a = calc_matrix_a(&params);

        // [rotation, shift, coupling] for flipped parameters with disorder
        let expected = [
            [0.0, 89.7, -4.24],
            [4.53778, 97.2, -2.82],
            [2.79248, 87.6, -2.27],
            [1.04718, 81.7, -3.82],
        ];

        let tol = 1e-1;

        for row in expected {
            let (rotation, expected_shift, expected_coupling) =
                (row[0], row[1], row[2]);

            let paf_coords = calc_paf_coords(rotation, &params, &matrix_a);
            let shift = calc_shift(&params, &paf_coords);
            let coupling = calc_dipolar_coupling(&params, &paf_coords);

            assert!(
                (shift - expected_shift).abs() < tol,
                "shift mismatch at rotation {}: got {}, expected {}",
                rotation,
                shift,
                expected_shift
            );

            assert!(
                (coupling - expected_coupling).abs() < tol,
                "coupling mismatch at rotation {}: got {}, expected {}",
                rotation,
                coupling,
                expected_coupling
            );
        }
    }
}
