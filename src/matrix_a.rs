use crate::params::Params;
use nalgebra::{Matrix3, Vector3};
use std::f64::consts::PI;

#[inline]
fn r1(angle: f64) -> Matrix3<f64> {
    let c = angle.cos();
    let s = angle.sin();
    Matrix3::new(1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c)
}

#[inline]
fn r3(angle: f64) -> Matrix3<f64> {
    let c = angle.cos();
    let s = angle.sin();
    Matrix3::new(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0)
}

/// Generate matrix A relating the principal axes frame (PAF) of the chemical shift
/// tensor to the helical axis frame (HAF).
///
/// # Arguments
/// * `params` - Fixed input parameters for the calculation
pub fn calc_matrix_a(params: &Params) -> Matrix3<f64> {
    // Precompute angles
    let e_alpha = PI - params.ca_c_n;
    let e_beta = PI - params.c_n_ca;
    let e_gamma = PI - params.n_ca_c;
    let omega = PI;

    // Compute v
    let v = Vector3::new(
        params.c_n * e_beta.cos()
            + params.ca_c * (e_alpha - e_beta).cos()
            + params.n_ca,
        params.c_n * -e_beta.sin() + params.ca_c * (e_alpha - e_beta).sin(),
        0.0,
    );

    // Combined rotation
    let c = r3(e_beta)
        * r1(omega)
        * r3(e_alpha)
        * r1(params.psi)
        * r3(e_gamma)
        * r1(params.phi);

    // Extract axis components
    let a = Vector3::new(
        c[(1, 2)] - c[(2, 1)],
        c[(2, 0)] - c[(0, 2)],
        c[(0, 1)] - c[(1, 0)],
    )
    .normalize();

    // Compute p
    let factor = 1.0 / (50.0_f64.to_radians().tan());
    let p = 0.5 * (factor * a.cross(&v) + v - a.dot(&v) * a);

    // Normalize r
    let r = -p.normalize();

    // Build HAF
    let haf = Matrix3::from_columns(&[r, a.cross(&r), a]).transpose();

    // Build R3 for final transform
    let r3 = r3(params.beta + params.ca_n_h);

    // Build permutation matrix for final transform
    let perm = Matrix3::new(0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    // Final matrix
    haf * r3 * perm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Params;

    fn mismatch(
        i: usize,
        j: usize,
        got: f64,
        want: f64,
        diff: f64,
        tol: f64,
    ) -> String {
        format!(
            "Mismatch at ({}, {}):\n\
            got = {:>10.6}\n\
            expected = {:>10.6}\n\
            diff = {:>10.6}\n\
            tol = {:>10.6}",
            i, j, got, want, diff, tol
        )
    }

    #[test]
    fn test_denny() {
        let params = Params::from_denny();
        let m = calc_matrix_a(&params);

        // Paper matrix (row-major)
        let paper = [
            [-0.83, 0.55, -0.09],
            [0.56, 0.80, -0.21],
            [-0.04, -0.22, -0.97],
        ];

        // Convert to nalgebra's column-major convention
        let expected = [
            [paper[0][0], paper[1][0], paper[2][0]],
            [paper[0][1], paper[1][1], paper[2][1]],
            [paper[0][2], paper[1][2], paper[2][2]],
        ];

        let tol = 0.01;

        for i in 0..3 {
            for j in 0..3 {
                let got = m[(i, j)];
                let want = expected[i][j];
                let diff = (got - want).abs();
                assert!(diff < tol, "{}", mismatch(i, j, got, want, diff, tol));
            }
        }
    }
}
