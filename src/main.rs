mod matrix_a;
mod paf;
mod params;
mod peak;
mod residue;
mod wave;

fn main() {
    let params = params::Params::new();
    let matrix_a = matrix_a::calc_matrix_a(&params);
    let rotation = 0.0; // Example rotation angle in radians
    let paf_coords = paf::calc_paf_coords(rotation, &params, &matrix_a);
    let shift = peak::calc_shift(&params, &paf_coords);
    let coupling = peak::calc_dipolar_coupling(&params, &paf_coords);
    let mut wave = wave::Wave::new(100); // Example with 100 points
    wave.update(&params, &matrix_a);

    /// SLN
    let mut params = params::Params::new();
    params.update_flip_angle(90.0);
    params.update_order_parameter(0.9);
    params.tilt = 24.6_f64.to_radians();
    params.rotation = 46.0_f64.to_radians();
    let mut residue = residue::Residue::new(-3, &params);
    residue.print_info();
    residue.calc_peak(&params, &matrix_a);
    residue.print_info();
    residue.set_exp_shift(135.499);
    residue.print_info();
    residue.set_exp_dipolar_coupling(1.348);
    residue.print_info();
    residue.unset_exp_shift();
    residue.print_info();

    println!("Matrix A: {:?}", matrix_a);
    println!("PAF Coordinates: {:?}", paf_coords);
    println!("Calculated Shift: {}", shift);
    println!("Calculated Dipolar Coupling: {}", coupling);
    println!("Wave Data: {:?}", wave.data);
}
