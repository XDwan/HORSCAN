use clap::Parser;
use anyhow::Result;

#[derive(Parser, Debug, Clone)]
#[command(
    name = "HORSCANv",
    version,
    // 简短描述 (显示在 -h 中)
    about = "Structure-aware hierarchical global aligner for centromeric regions.",
    // 详细描述 (显示在 --help 中)
    long_about = "HORSCANv: A Two-Level Hierarchical Aligner (v2.2)

Performs global alignment on genomic sequences with Higher-Order Repeat (HOR) structures.
It uses a 'State-Aware' algorithm that combines:
  1. Micro-alignment: Affine Gap Needleman-Wunsch (Gotoh) for monomers.
  2. Macro-alignment: Structural jumps between HOR units based on pre-computed anchors.

NOTE: All scoring inputs should be POSITIVE magnitudes. Penalties (mismatch, gaps) are automatically negated internally."
)]
pub struct Args {
    // --- Input Files ---
    #[arg(short = 's', long, help_heading = "Input/Output")]
    /// Path to the Source monomer BED file
    pub source: String,

    #[arg(short = 't', long, help_heading = "Input/Output")]
    /// Path to the Target monomer BED file
    pub target: String,

    #[arg(long = "source-hor", help_heading = "Input/Output")]
    /// Path to the Source HOR annotation BED file
    pub source_hor: String,

    #[arg(long = "target-hor", help_heading = "Input/Output")]
    /// Path to the Target HOR annotation BED file
    pub target_hor: String,

    #[arg(short = 'o', long, default_value = "out", help_heading = "Input/Output")]
    /// Prefix for output files (e.g., "results/chr1")
    pub output_prefix: String,

    // --- Scoring Parameters ---
    
    /// Monomer alignment scoring parameters (Affine only).
    ///
    /// USAGE:
    ///   -m <MATCH> <MISMATCH> <OPEN> <EXTEND>
    ///
    /// EXAMPLE:
    ///   -m 10 4 2 2  (Match=+10, Mismatch=-4, Gap=-2-2*k)
    #[arg(short = 'm', num_args = 4, required = true, help_heading = "Scoring")]
    pub mon_scores: Vec<i32>,

    /// HOR structural scoring parameters (Affine + label mismatch penalty).
    ///
    /// USAGE:
    ///   -h <MATCH> <MISMATCH> <OPEN> <EXTEND>
    ///
    /// EXAMPLE:
    ///   -h 20 10 1 1  (Match=+20, Mismatch=-10, Gap=-1-1*k)
    #[arg(short = 'h', long = "hor", num_args = 4, required = true, help_heading = "Scoring")]
    pub hor_scores: Vec<i32>,

    // --- Heuristics ---

    /// Band width for pre-filtering HOR pairs (0 to disable).
    /// Optimization: Only HOR pairs within this diagonal distance are considered for macro-jumps.
    #[arg(long = "hor-pair-band", default_value_t = 0, help_heading = "Heuristics")]
    pub hor_pair_band: usize,

    /// Allow relaxed HOR pairing (label equivalence) instead of exact match.
    #[arg(long, default_value_t = false, help_heading = "Heuristics")]
    pub relaxed: bool,

    // --- Misc ---

    /// Emit detailed debug information (e.g., .tsv logs).
    #[arg(long, default_value_t = false, help_heading = "Misc")]
    pub debug: bool,
}

impl Args {
    pub fn parse() -> Self {
        <Self as Parser>::parse()
    }
}

pub fn validate_args(a: &Args) -> Result<()> {
    // 1. Validate Scoring Mode (Affine only)
    anyhow::ensure!(
        a.mon_scores.len() == 4,
        "Invalid number of scores for -m. Expected 4 (Affine: MATCH MISMATCH OPEN EXTEND)."
    );
    anyhow::ensure!(
        a.hor_scores.len() == 4,
        "Invalid number of scores for -h. Expected 4 (Affine: MATCH MISMATCH OPEN EXTEND)."
    );

    // 2. Validate Positive Magnitudes
    anyhow::ensure!(
        a.mon_scores.iter().all(|&x| x >= 0), 
        "All micro scores must be non-negative integers (magnitudes)."
    );
    anyhow::ensure!(
        a.hor_scores.iter().all(|&x| x >= 0), 
        "All macro scores must be non-negative integers (magnitudes)."
    );
    
    // 3. Logic Check: HOR extension shouldn't be more expensive than Monomer extension
    // If HOR extend is more expensive, the algorithm might never choose the HOR path for deletions.
    let mon_gap_extend = a.mon_scores[3];
    let hor_gap_extend = a.hor_scores[3];
    
    anyhow::ensure!(
        hor_gap_extend <= mon_gap_extend, 
        "Logic Error: HOR gap extension cost ({}) cannot be higher than Monomer gap extension cost ({}).",
        hor_gap_extend, mon_gap_extend
    );

    Ok(())
}