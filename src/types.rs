#![allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ModeMon {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

#[derive(Debug, Clone)]
pub struct ModeHor {
    pub mon: ModeMon,
    pub relaxed: bool,
    // Macro level label scores
    pub hor_match_score: i32,      // +S_match^{macro}
    pub hor_mismatch_score: i32,   // -P_mismatch^{macro}
    // Macro level affine gap for HOR indels
    pub hor_gap_open: i32,         // -P_open^{macro}
    pub hor_gap_extend: i32,       // -P_extend^{macro}
}

#[derive(Debug, Clone)]
pub struct MonRow {
    pub sample: String,
    pub start: i64,
    pub end: i64,
    pub mon: String,
}

#[derive(Debug, Clone)]
pub struct HorRow {
    pub sample: String,
    pub start: i64,
    pub end: i64,
    pub hor_label: String,
}

#[derive(Debug, Clone)]
pub struct HorUnit {
    pub label: String,
    pub start_mi: usize, // 0-based index relative to original monomer list WITHOUT sentinel (sentinel lives at global index 0)
    pub end_mi: usize,   // inclusive (same coordinate system as start_mi)
}
impl HorUnit { pub fn len(&self) -> usize { self.end_mi - self.start_mi + 1 } }

#[derive(Debug, Clone)]
pub struct MonAlignmentRow {
    pub src_sample: String,
    pub src_start: i64,
    pub src_end: i64,
    pub src_mon: String,
    pub tgt_sample: String,
    pub tgt_start: i64,
    pub tgt_end: i64,
    pub tgt_mon: String,
    // align_type: "MTH" | "MIS" | "INS" | "DEL"
    pub align_type: String,
    // Added: original indices (including sentinel=0). 0 means gap on that side.
    pub src_mi: usize,
    pub tgt_mi: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum StepCode { Diag=0, Up=1, Left=2, HorMatch=3, HorDel=4, HorIns=5 }

#[derive(Debug, Clone)]
pub struct DebugStep { pub i: usize, pub j: usize, pub step: StepCode }

#[derive(Debug, Clone)]
pub struct AlignOutputs {
    pub rows: Vec<MonAlignmentRow>,
    pub debug: Vec<DebugStep>,
    pub hor_events: Vec<HorEvent>,
}

// Tagging: map each monomer index -> Option<HOR id>
#[derive(Debug, Clone)]
pub struct HorTagging {
    pub mon_to_hor: Vec<Option<usize>>, // length = mon_count+1 (index 0 sentinel placeholder)
    pub units: Vec<HorUnit>,
}

#[derive(Debug, Clone)]
pub struct HorEvent {
    pub event_type: String, // MATCH | INS | DEL (INS = source insertion vs target gap; DEL = source deletion vs target present)
    // source HOR info (if present)
    pub src_label: String,
    pub src_hor_id: Option<usize>,
    pub src_start_mi: usize,
    pub src_end_mi: usize,
    pub src_start_bp: i64,
    pub src_end_bp: i64,
    // target HOR info (if present)
    pub tgt_label: String,
    pub tgt_hor_id: Option<usize>,
    pub tgt_start_mi: usize,
    pub tgt_end_mi: usize,
    pub tgt_start_bp: i64,
    pub tgt_end_bp: i64,
    // scoring summary
    pub score: i32,
    // micro alignment stats within this HOR event
    pub mth: usize,
    pub mis: usize,
    pub ins: usize,
    pub del: usize,
    // shift statistics (only meaningful for MATCH events with >=1 aligned monomer pair)
    pub mean_shift: Option<f64>,
    pub shift_std: Option<f64>,
    pub shift_consistent: bool,
}

// Removed auxiliary output structs: ShiftRow, UnalignedRow, StatsSummary to slim data model.
