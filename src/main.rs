mod params;
mod matrix_a;
use params::Params;

fn main() {
    let params = Params::new();
    let matrix_a = matrix_a::generate_matrix_a(&params);
    println!("Matrix A: {:?}", matrix_a);
}
