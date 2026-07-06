mod matrix_a;
mod paf;
mod params;
mod peak;
mod residue;
mod wave;
mod wheel;

use residue::ResidueData;
use wheel::Wheel;

fn main() {

    /// Create a wheel from SLN experimental data.
    /// https://github.com/weberdak/pisa.py/blob/master/examples/sarcolipin/sln_explore_log.dat
    fn create_sln_wheel() -> wheel::Wheel {
        let mut params = params::Params::new();
        params.update_flip_angle(90.0);
        params.update_order_parameter(0.9);
        params.tilt = 24.6_f64.to_radians();
        params.rotation = 46.0_f64.to_radians();

        let data = vec![
            ResidueData::new(-5, Some([b'M']), None, None, true), // M1
            ResidueData::new(-4, Some([b'G']), None, None, true), // G2
            ResidueData::new(-3, Some([b'I']), Some(135.499), Some(1.348), true), // I3
            ResidueData::new(-2, Some([b'N']), None, None, true), // N4
            ResidueData::new(-1, Some([b'T']), None, None, true), // T5
            ResidueData::new(0, Some([b'R']), None, None, true),  // R6
            ResidueData::new(1, Some([b'E']), None, None, true),  // E7
            ResidueData::new(2, Some([b'L']), Some(94.699), Some(1.011), true), // L8
            ResidueData::new(3, Some([b'F']), Some(80.834), Some(2.624), true), // F9
            ResidueData::new(4, Some([b'L']), Some(101.172), Some(4.048), true), // L10
            ResidueData::new(5, Some([b'N']), Some(108.467), Some(2.429), true), // N11
            ResidueData::new(6, Some([b'F']), Some(86.527), Some(1.813), true), // F12
            ResidueData::new(7, Some([b'T']), Some(77.984), Some(3.943), true), // T13
            ResidueData::new(8, Some([b'I']), Some(101.172), Some(4.048), false), // I14
            ResidueData::new(9, Some([b'V']), Some(97.199), Some(2.116), false), // V15
            ResidueData::new(10, Some([b'L']), Some(78.544), Some(2.753), false), // L16
            ResidueData::new(11, Some([b'I']), Some(83.056), Some(4.483), false), // I17
            ResidueData::new(12, Some([b'T']), Some(99.822), Some(3.357), false), // T18
            ResidueData::new(13, Some([b'V']), Some(89.691), Some(1.953), false), // V19
            ResidueData::new(14, Some([b'I']), Some(77.359), Some(3.494), false), // I20
            ResidueData::new(15, Some([b'L']), Some(93.226), Some(4.536), false), // L21
            ResidueData::new(16, Some([b'M']), Some(98.304), Some(2.661), false), // M22
            ResidueData::new(17, Some([b'W']), Some(85.096), Some(2.211), false), // W23
            ResidueData::new(18, Some([b'L']), Some(80.493), Some(3.964), false), // L24
            ResidueData::new(19, Some([b'L']), Some(101.172), Some(4.048), false), // L25
            ResidueData::new(20, Some([b'V']), Some(94.381), Some(2.116), false), // V26
            ResidueData::new(21, Some([b'R']), Some(78.544), Some(2.753), false), // R27
            ResidueData::new(22, Some([b'S']), Some(83.277), Some(4.144), true), // S28
            ResidueData::new(23, Some([b'Y']), Some(114.132), Some(1.075), true), // Y29
            ResidueData::new(24, Some([b'Q']), None, None, true), // Q30
            ResidueData::new(25, Some([b'Y']), Some(110.256), Some(1.606), true), // Y31
        ];
        Wheel::from_data(params, 6, 1, data)
    }

    // Print Experimental and Simulated peaks for the SLN wheel
    //let wheel = create_sln_wheel();
    //println!("{}", wheel.as_csv(3));
    let params = params::Params::new();
    let mut wheel = Wheel::new(params, 5, 1, 1);
    wheel.update_exp_shift(1, Some(205.0));
    wheel.update_exp_dipolar_coupling(1, Some(10.0));
    println!("{}", wheel.as_csv(6));
    println!("{}", wheel.rmsd);
    wheel.update_coupling_constant(15.0);
    println!("{}", wheel.as_csv(6));
    println!("{}", wheel.rmsd);
}
