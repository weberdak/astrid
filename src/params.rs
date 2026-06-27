/// Parameter class for all user-defines inputs.
pub struct Params {
    /// Helical tilt angle in radians (0-pi/2)
    pub tilt: f64,
    /// Helical rotation angle at reference residue in radians (0-pi/2)
    pub rotation: f64,
    /// Order parameter of the helix (0-1)
    pub order: f64,
    /// Flip angle of the helix in radians (0 or pi/2)
    pub flip: f64,
    /// Phi dihedral angle for an ideal helix in radians (-pi to pi)
    pub phi: f64,
    /// Psi dihedral angle for an ideal helix in radians (-pi to pi)
    pub psi: f64,
    /// Beta rotation of the 15N PAS with respect to the amide plane in radians (0-pi)
    pub beta: f64,
    /// CA-C-N bond angle in radians
    pub ca_c_n: f64,
    /// C-N-CA bond angle in radians
    pub c_n_ca: f64,
    /// N-CA-C bond angle in radians
    pub n_ca_c: f64,
    /// CA-N-H bond angle in radians
    pub ca_n_h: f64,
    /// CA-C bond length in angstroms
    pub ca_c: f64,
    /// C-N bond length in angstroms
    pub c_n: f64,
    /// N-CA bond length in angstroms
    pub n_ca: f64,
}

impl Params {
    /// Create with default parameters from Weber et al. (2020)
    /// https://doi.org/10.1093/bioinformatics/btaa019
    pub fn new() -> Self {
        Self { 
            tilt: 20.0_f64.to_radians(),
            rotation: 0.0_f64.to_radians(),
            order: 1.0,
            flip: 0.0_f64.to_radians(),
            phi: -63.0_f64.to_radians(),
            psi: -42.0_f64.to_radians(),
            beta: 17.0_f64.to_radians(),
            ca_c_n: 117.5_f64.to_radians(),
            c_n_ca: 124.0_f64.to_radians(),
            n_ca_c: 107.4_f64.to_radians(),
            ca_n_h: 116.0_f64.to_radians(),
            ca_c: 1.52,
            c_n: 1.35,
            n_ca: 1.45,
        }
    }

    /// Create with parameters from original derivation by Denny et al. (2001)
    /// https://doi.org/10.1006/jmre.2001.2405
    pub fn from_denny() -> Self {
        Self {
            tilt: 35.0_f64.to_radians(),
            rotation: -10.0_f64.to_radians(),
            order: 1.0,
            flip: 0.0_f64.to_radians(),
            phi: -65.0_f64.to_radians(),
            psi: -40.0_f64.to_radians(),
            beta: 17.0_f64.to_radians(),
            ca_c_n: 115.0_f64.to_radians(),
            c_n_ca: 121.0_f64.to_radians(),
            n_ca_c: 110.0_f64.to_radians(),
            ca_n_h: 117.0_f64.to_radians(),
            ca_c: 1.53,
            c_n: 1.34,
            n_ca: 1.45,
        }
    }
}
