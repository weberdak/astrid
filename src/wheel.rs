





pub struct Wheel {
    /// Number of residues in the sequence
    pub num_residues: usize,
    /// Simulated peak for each residue in the sequence, (shift, coupling)
    pub sim_data: Vec<(f64, f64)>,
    /// Experimental peak for each residue in the sequence, (shift, coupling)
    pub exp_data: Vec<(f64, f64)>,
    /// RMSD between simulated and experimental peaks
    pub rmsd: f64,
    /// Sequence of residues (e.g., "ACDEFGHIKLMNPQRSTVWY")
    pub sequence: String,
    /// Start index of the residue in the sequence for the simulation
    pub start_index: usize,
    /// Index of reference residue in the sequence for the simulation
    pub reference_index: usize,
}