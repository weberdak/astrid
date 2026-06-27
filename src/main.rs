mod helix;
mod matrix_a;
use helix::Helix;

fn main() {
    let helix = Helix::new();
    println!("Matrix A: {:?}", helix.matrix_a);
}
