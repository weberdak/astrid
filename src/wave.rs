use crate::paf::calc_paf_coords;
use crate::params::Params;
use crate::peak::{calc_dipolar_coupling, calc_shift};
use nalgebra::Matrix3;
use std::f64::consts::PI;

/// Stores a continuous wave of calculated shifts and couplings for a given number
/// of points in one full rotation (0 to 2π radians).
pub struct Wave {
    pub data: Vec<(f64, f64, f64)>, // (rotation, shift, coupling)
    pub num_points: usize,
}

impl Wave {
    pub fn new(num_points: usize) -> Self {
        Wave {
            data: Vec::with_capacity(num_points),
            num_points,
        }
    }

    pub fn update(&mut self, params: &Params, matrix_a: &Matrix3<f64>) {
        let step = 2.0 * PI / self.num_points as f64;

        self.data = (0..self.num_points)
            .map(|i| {
                let rotation = i as f64 * step;
                let paf = calc_paf_coords(rotation, params, matrix_a);
                let shift = calc_shift(params, &paf);
                let coupling = calc_dipolar_coupling(params, &paf);
                (rotation, shift, coupling)
            })
            .collect();
    }
}
