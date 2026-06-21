mod matrix_a;
mod helix;
use helix::Helix;

fn main() {
    let helix = Helix::new();
    println!("Matrix A: {:?}", helix.matrix_a);
}
