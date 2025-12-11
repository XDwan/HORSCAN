use crate::types::*;

#[inline]
pub fn omega_monomer(a: &str, b: &str, mode: &ModeMon) -> i32 {
    if a == b { mode.match_score } else { mode.mismatch_score }
}

#[inline]
pub fn hor_gap(len_mons: usize, mode: &ModeHor) -> i32 {
    // HORGap(H) = open + len * extend  (macro affine gap)
    mode.hor_gap_open + (len_mons as i32) * mode.hor_gap_extend
}
