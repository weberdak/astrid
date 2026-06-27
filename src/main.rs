mod matrix_a;
mod paf;
mod params;
mod peak;
use params::Params;

fn main() {
    let params = Params::new();
    let matrix_a = matrix_a::calc_matrix_a(&params);
    let rotation = 0.0; // Example rotation angle in radians
    let paf_coords = paf::calc_paf_coords(rotation, &params, &matrix_a);
    let shift = peak::calc_shift(&params, &paf_coords);
    let coupling = peak::calc_dipolar_coupling(&params, &paf_coords);
    println!("Matrix A: {:?}", matrix_a);
    println!("PAF Coordinates: {:?}", paf_coords);
    println!("Calculated Shift: {}", shift);
    println!("Calculated Dipolar Coupling: {}", coupling);
}
