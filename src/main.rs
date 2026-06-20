mod matrix_a;
mod helix;
use helix::Helix;

fn main() {
    let helix = Helix::new();
    println!("{:?}", helix.matrix_a);
}
