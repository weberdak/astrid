
/// Fixed input parameters for calculating a PISA wheel.
pub struct Params {
    /// Helical tilt angle in degrees (0-90)
    pub tau_deg: f64,
    /// Helical rotation angle in radians (0-pi/2)
    pub rho_rad: f64,
    /// Azimuthal angle of reference residue in degrees (0-360)
    pub rho0_deg: f64,
    /// Order parameter of the helix (0-1)
    pub order: f64,
    /// Flip angle of the helix in degrees (0-90)
    pub flip_deg: f64,
    /// Phi dihedral angle for an ideal helix in degrees (-180 to 180)
    pub phi_deg: f64,
    /// Psi dihedral angle for an ideal helix in degrees (-180 to 180)
    pub psi_deg: f64,
    /// Beta rotation of the 15N PAS with respect to the amide plane in degrees
    pub beta_deg: f64,
    /// Maximum 15N-1H dipolar coupling in kHz
    pub dcMax: f64,
    /// CA,C,N bond angle (radians)
    pub aCaCN: f64,

}
