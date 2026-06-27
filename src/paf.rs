use crate::params::Params;
use nalgebra::{Matrix3, Vector3};

/// Calculate coordinates of B0 in Principal Axes Frame (PAF) of a residue at a given
/// azimuthal rotation angle along the helix.
///
/// # Arguments
/// * `rotation` - Azimuthal angle of residue in radians
/// * `params` - The parameters for the calculation
/// * `matrix_a` - Rotation matrix relating the principal axes frame (PAF) of the
///   chemical shift tensor to the helical axis frame (HAF)
pub fn calc_paf_coords(
    rotation: f64,
    params: &Params,
    matrix_a: &Matrix3<f64>,
) -> Vector3<f64> {
    let matrix_x = Vector3::new(
        rotation.cos() * params.tilt.sin(),
        rotation.sin() * params.tilt.sin(),
        params.tilt.cos(),
    );
    matrix_a.transpose() * matrix_x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix_a::calc_matrix_a;
    use crate::params::Params;

    #[test]
    fn test_calc_paf_coords() {
        let params = Params::new();
        let matrix_a = calc_matrix_a(&params);

        // Expect [-0.37897365, 0.02097897, -0.92516964] @ 0.0 rotation
        let paf_coords = calc_paf_coords(0.0, &params, &matrix_a);
        let expected = Vector3::new(-0.37897365, 0.02097897, -0.92516964);

        let tol = 1e-4;
        assert!((paf_coords.x - expected.x).abs() < tol);
        assert!((paf_coords.y - expected.y).abs() < tol);
        assert!((paf_coords.z - expected.z).abs() < tol);
    }
}
