/// Parameter class for all fixed user-defined inputs.
#[derive(Clone, Debug)]
pub struct Params {
    /// Helical tilt angle in radians (0-π/2)
    pub tilt: f64,
    /// Helical rotation angle at reference residue in radians (0-π/2)
    pub rotation: f64,
    /// Flip angle of bicelle in radians (0-π/2)
    pub flip: f64,
    /// Order parameter of the helix (0-1)
    pub order: f64,
    /// Scaling factor to account for bicelle flip angle and helical order parameter
    pub scalar: f64,
    /// Phi dihedral angle for an ideal helix in radians (-π to π)
    pub phi: f64,
    /// Psi dihedral angle for an ideal helix in radians (-π to π)
    pub psi: f64,
    /// Beta rotation of the 15N PAS with respect to the amide plane in radians (0-π)
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
    /// Smallest component of the 15N CSA tensor in ppm
    pub delta_xx: f64,
    /// Largest component of the 15N CSA tensor in ppm
    pub delta_yy: f64,
    /// Intermediate component of the 15N CSA tensor in ppm
    pub delta_zz: f64,
    /// Isotropic chemical shift calculated from the 15N CSA tensor in ppm
    pub shift_iso: f64,
    /// Dipolar coupling constant for the 15N-1H bond in kHz
    pub coupling: f64,
}

impl Params {
    /// Create with default parameters from Weber et al. (2020)
    /// https://doi.org/10.1093/bioinformatics/btaa019
    pub fn new() -> Self {
        let flip = 0.0_f64.to_radians();
        let order = 1.0;
        let delta_xx = 57.3;
        let delta_yy = 81.2;
        let delta_zz = 228.1;
        Self {
            tilt: 20.0_f64.to_radians(),
            rotation: 0.0_f64.to_radians(),
            flip,
            order,
            scalar: 0.5 * (3.0 * flip.cos().powi(2) - 1.0) * order,
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
            delta_xx,
            delta_yy,
            delta_zz,
            shift_iso: (delta_xx + delta_yy + delta_zz) / 3.0,
            coupling: 10.735,
        }
    }

    /// Create with parameters from original derivation by Denny et al. (2001)
    /// https://doi.org/10.1006/jmre.2001.2405
    pub fn from_denny() -> Self {
        let flip = 0.0_f64.to_radians();
        let order = 1.0;
        let delta_xx = 31.0;
        let delta_yy = 55.0;
        let delta_zz = 202.0;
        Self {
            tilt: 35.0_f64.to_radians(),
            rotation: -10.0_f64.to_radians(),
            flip,
            order,
            scalar: 0.5 * (3.0 * flip.cos().powi(2) - 1.0) * order,
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
            delta_xx,
            delta_yy,
            delta_zz,
            shift_iso: (delta_xx + delta_yy + delta_zz) / 3.0,
            coupling: 11.335,
        }
    }

    /// Set flip angle in degrees and recalculate the scalar factor
    pub fn update_flip_angle(&mut self, flip_angle_deg: f64) {
        self.flip = flip_angle_deg.to_radians();
        self.scalar = 0.5 * (3.0 * self.flip.cos().powi(2) - 1.0) * self.order;
    }

    /// Set order parameter and recalculate the scalar factor
    pub fn update_order_parameter(&mut self, order: f64) {
        self.order = order;
        self.scalar = 0.5 * (3.0 * self.flip.cos().powi(2) - 1.0) * self.order;
    }

    /// Set the δxx component of the 15N CSA tensor and recalculate the isotropic shift
    pub fn update_delta_xx(&mut self, delta_xx: f64) {
        self.delta_xx = delta_xx;
        self.shift_iso = (self.delta_xx + self.delta_yy + self.delta_zz) / 3.0;
    }

    /// Set the δyy component of the 15N CSA tensor and recalculate the isotropic shift
    pub fn update_delta_yy(&mut self, delta_yy: f64) {
        self.delta_yy = delta_yy;
        self.shift_iso = (self.delta_xx + self.delta_yy + self.delta_zz) / 3.0;
    }

    /// Set the δzz component of the 15N CSA tensor and recalculate the isotropic shift
    pub fn update_delta_zz(&mut self, delta_zz: f64) {
        self.delta_zz = delta_zz;
        self.shift_iso = (self.delta_xx + self.delta_yy + self.delta_zz) / 3.0;
    }
}
