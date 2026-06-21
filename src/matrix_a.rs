use crate::helix::Helix;
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
pub fn generate_matrix_a(h: &Helix) -> Matrix3<f64> {
    // Precompute angles
    let e_alpha = PI - h.ca_c_n;
    let e_beta = PI - h.c_n_ca;
    let e_gamma = PI - h.n_ca_c;
    let omega = PI;

    // Compute v
    let v = Vector3::new(
        h.c_n * e_beta.cos() + h.ca_c * (e_alpha - e_beta).cos() + h.n_ca,
        h.c_n * -e_beta.sin() + h.ca_c * (e_alpha - e_beta).sin(),
        0.0,
    );

    // Combined rotation
    let c = r3(e_beta) * r1(omega) * r3(e_alpha) * r1(h.psi) * r3(e_gamma) * r1(h.phi);

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
    let r3 = r3(h.beta + h.ca_n_h);

    // Build permutation matrix for final transform
    let perm = Matrix3::new(0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    // Final matrix
    haf * r3 * perm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::helix::Helix;

    fn mismatch(i: usize, j: usize, got: f64, want: f64, diff: f64, tol: f64) -> String {
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
        let helix = Helix::from_denny();
        let m = generate_matrix_a(&helix);

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
