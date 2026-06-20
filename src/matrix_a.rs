use crate::helix::Helix;
use nalgebra::{Matrix3, Vector3};
use std::f64::consts::PI;

#[inline]
fn r1(angle: f64) -> Matrix3<f64> {
    let c = angle.cos();
    let s = angle.sin();
    Matrix3::new(1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c)
}

#[inline]
fn r3(angle: f64) -> Matrix3<f64> {
    let c = angle.cos();
    let s = angle.sin();
    Matrix3::new(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0)
}

/// Generate the transformation matrix for a given helix.
pub fn generate_matrix_a(h: &Helix) -> Matrix3<f64> {
    // Precompute angles
    let e_alpha = PI - h.ca_c_n;
    let e_beta = PI - h.c_n_ca;
    let e_gamma = PI - h.n_ca_c;
    let omega = 180.0_f64.to_radians();

    // Compute v
    let v = Vector3::new(
        h.c_n * e_beta.cos() + h.ca_c * (e_alpha - e_beta).cos() + h.n_ca,
        h.c_n * -e_beta.sin() + h.ca_c * (e_alpha - e_beta).sin(),
        0.0,
    );

    // Rotation matrices
    let r1_phi = r1(h.phi);
    let r3_gamma = r3(e_gamma);
    let r1_psi = r1(h.psi);
    let r3_alpha = r3(e_alpha);
    let r1_omega = r1(omega);
    let r3_beta = r3(e_beta);

    // Combined rotation
    let c = r1_phi * r3_gamma * r1_psi * r3_alpha * r1_omega * r3_beta;

    // Extract axis components
    let a = Vector3::new(
        c[(2, 1)] - c[(1, 2)],
        c[(0, 2)] - c[(2, 0)],
        c[(1, 0)] - c[(0, 1)],
    )
    .normalize();

    // Compute p
    let factor = 1.0 / (50.0_f64.to_radians().tan());
    let p = 0.5 * (factor * a.cross(&v) + v - a.dot(&v) * a);

    // Normalize r
    let r = -p.normalize();

    // Build HAF
    let haf = Matrix3::from_columns(&[r, a.cross(&r), a]);

    // Build R3 for final transform
    let angle = h.beta + h.ca_n_h;
    let r3 = Matrix3::new(
        angle.cos(),
        -angle.sin(),
        0.0,
        angle.sin(),
        angle.cos(),
        0.0,
        0.0,
        0.0,
        1.0,
    );

    // Build permutation matrix for final transform
    let perm = Matrix3::new(0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);

    // Final matrix
    perm * r3 * haf
}
