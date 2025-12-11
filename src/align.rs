use crate::{errors::HorScanError, score::*, types::*};
use std::collections::HashMap;

fn idx(i: usize, j: usize, n: usize) -> usize {
    i * (n + 1) + j
}

/// Trace codes for hierarchical affine DP (only for M state)
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TraceCode {
    Empty   = 0, // unused / boundary
    Diag    = 1, // micro from M
    FromIx  = 2, // micro from Ix
    FromIy  = 3, // micro from Iy
    HorJump = 4, // any HOR macro-jump (MATCH / INS / DEL)
}

// Standard global NW for internal HOR alignment (score only + backtrace)
#[derive(Clone)]
pub struct MiniAlign {
    pub score: i32,
    pub path: Vec<(Option<usize>, Option<usize>)>,
}

/// Full Needleman-Wunsch producing score + path (O(mn) space).
/// Updated to use Affine Gap Penalty (Gotoh's algorithm) with lazy backtrace.
pub fn nw_global_mini(src: &[String], tgt: &[String], mode: &ModeMon) -> MiniAlign {
    let m = src.len();
    let n = tgt.len();
    let neg_inf = i32::MIN / 4;
    
    // We need 3 matrices: M, Ix, Iy
    // M: Best score ending in alignment (Match/Mismatch)
    // Ix: Best score ending in insertion in Target (vertical gap)
    // Iy: Best score ending in insertion in Source (horizontal gap)
    
    let mut mat_m = vec![neg_inf; (m + 1) * (n + 1)];
    let mut mat_ix = vec![neg_inf; (m + 1) * (n + 1)];
    let mut mat_iy = vec![neg_inf; (m + 1) * (n + 1)];
    
    // Initialization
    mat_m[idx(0, 0, n)] = 0;
    // Ix/Iy at (0,0) remain neg_inf (invalid to end in gap extension at start)
    
    // First col (gap in target)
    for i in 1..=m {
        // M[i][0] = Gap(i) = open + i*extend
        mat_m[idx(i, 0, n)] = mode.gap_open + (i as i32) * mode.gap_extend;
        // Ix[i][0]: valid to be in gap state.
        // If we define Ix as "end in gap extension", then Ix[i][0] should ideally equal M[i][0].
        // Let's set it consistent with recurrence:
        // Ix[i][0] = max(M[i-1][0]+open+ext, Ix[i-1][0]+ext). 
        // For i=1: M[0][0]+open+ext = open+ext.
        // NOTE: We initialize Ix to the gap score so that subsequent Ix calculations can extend from it.
        // If we set it to neg_inf, we would pay the open penalty at every step along the boundary.
        mat_ix[idx(i, 0, n)] = mat_m[idx(i, 0, n)];
        // Iy[i][0] = neg_inf
    }
    
    // First row (gap in source)
    for j in 1..=n {
        mat_m[idx(0, j, n)] = mode.gap_open + (j as i32) * mode.gap_extend;
        mat_iy[idx(0, j, n)] = mat_m[idx(0, j, n)];
        // Ix[0][j] = neg_inf
    }
    
    for i in 1..=m {
        for j in 1..=n {
            let idx_curr = idx(i, j, n);
            let idx_diag = idx(i-1, j-1, n);
            let idx_up   = idx(i-1, j, n);
            let idx_left = idx(i, j-1, n);
            
            // 1. Update Ix (Vertical)
            let m_up = mat_m[idx_up];
            let ix_up = mat_ix[idx_up];
            let score_ix = (m_up + mode.gap_open + mode.gap_extend).max(ix_up + mode.gap_extend);
            mat_ix[idx_curr] = score_ix;
            
            // 2. Update Iy (Horizontal)
            let m_left = mat_m[idx_left];
            let iy_left = mat_iy[idx_left];
            let score_iy = (m_left + mode.gap_open + mode.gap_extend).max(iy_left + mode.gap_extend);
            mat_iy[idx_curr] = score_iy;
            
            // 3. Update M (Match)
            let score_base = omega_monomer(&src[i-1], &tgt[j-1], mode);
            let best_diag = mat_m[idx_diag].max(mat_ix[idx_diag]).max(mat_iy[idx_diag]);
            mat_m[idx_curr] = best_diag + score_base;
        }
    }
    
    let final_score = mat_m[idx(m, n, n)].max(mat_ix[idx(m, n, n)]).max(mat_iy[idx(m, n, n)]);
    
    // Lazy Backtrace
    let mut path = Vec::new();
    let (mut i, mut j) = (m, n);
    
    // Determine start state
    let mut state = if i == 0 && j > 0 { 2 } // Boundary: must be Iy
               else if j == 0 && i > 0 { 1 } // Boundary: must be Ix
               else if mat_m[idx(m, n, n)] == final_score { 0 } // M
               else if mat_ix[idx(m, n, n)] == final_score { 1 } // Ix
               else { 2 }; // Iy
               
    while i > 0 || j > 0 {
        match state {
            0 => { // M state: came from (i-1, j-1)
                // Boundary guard: avoid using i-1 or j-1 when at edges
                if i == 0 && j > 0 { state = 2; continue; }
                if j == 0 && i > 0 { state = 1; continue; }
                if i == 0 && j == 0 { break; }
                path.push((Some(i - 1), Some(j - 1)));
                let idx_diag = idx(i-1, j-1, n);
                let score_needed = mat_m[idx(i, j, n)] - omega_monomer(&src[i-1], &tgt[j-1], mode);
                
                i -= 1;
                j -= 1;

                if i == 0 && j > 0 { state = 2; }
                else if j == 0 && i > 0 { state = 1; }
                else if i == 0 && j == 0 { state = 0; } // Done
                else if mat_m[idx_diag] == score_needed { state = 0; }
                else if mat_ix[idx_diag] == score_needed { state = 1; }
                else { state = 2; } // assume iy matches if others don't (tie-break)
            },
            1 => { // Ix state: came from (i-1, j)
                // Boundary guard for top row
                if i == 0 {
                    if j > 0 { state = 2; continue; } else { break; }
                }
                path.push((Some(i - 1), None));
                let idx_up = idx(i-1, j, n);
                let curr_score = mat_ix[idx(i, j, n)];
                // Did we open or extend?
                // Open: M[i-1][j] + open + extend
                // Extend: Ix[i-1][j] + extend
                let score_ext = mat_ix[idx_up] + mode.gap_extend;
                if score_ext == curr_score { state = 1; }
                else { state = 0; } // Must be open from M
                i -= 1;
            },
            2 => { // Iy state: came from (i, j-1)
                // Boundary guard for left column
                if j == 0 {
                    if i > 0 { state = 1; continue; } else { break; }
                }
                path.push((None, Some(j - 1)));
                let idx_left = idx(i, j-1, n);
                let curr_score = mat_iy[idx(i, j, n)];
                // Extend: Iy[i][j-1] + extend
                let score_ext = mat_iy[idx_left] + mode.gap_extend;
                if score_ext == curr_score { state = 2; }
                else { state = 0; }
                j -= 1;
            },
            _ => unreachable!()
        }
    }
    path.reverse();
    MiniAlign {
        score: final_score,
        path,
    }
}

/// Linear-space (two-row) global NW computing ONLY the final score.
/// Updated to use Affine Gap Penalty (Gotoh's algorithm).
pub fn nw_global_mini_score(src: &[String], tgt: &[String], mode: &ModeMon) -> i32 {
    let m = src.len();
    let n = tgt.len();
    // Special cases for empty sequences using affine logic
    // Gap(k) = open + k*extend
    if m == 0 { 
        return if n == 0 { 0 } else { mode.gap_open + (n as i32) * mode.gap_extend };
    }
    if n == 0 { 
        return if m == 0 { 0 } else { mode.gap_open + (m as i32) * mode.gap_extend };
    }

    let neg_inf = i32::MIN / 4;

    // We need to store the previous row's M and Ix scores
    // prev_M[j] stores M[i-1][j]
    // prev_Ix[j] stores Ix[i-1][j]
    let mut prev_m = vec![neg_inf; n + 1];
    let mut prev_ix = vec![neg_inf; n + 1];
    
    // Initialize 0-th row (i=0)
    prev_m[0] = 0;
    prev_ix[0] = neg_inf; // Cannot start with insertion relative to Target at (0,0)
    // Iy at (0,0) is also neg_inf implicitly

    // Fill 0-th row boundaries (gap in Target, i.e., Source=0, Target=j)
    // This corresponds to Iy state theoretically, but M tracks the valid path
    for j in 1..=n {
        prev_m[j] = mode.gap_open + (j as i32) * mode.gap_extend;
        // prev_ix[j] remains neg_inf as we can't have vertical gap at row 0
    }

    for i in 1..=m {
        // Start of row i
        // curr_iy tracks Iy[i][j] (gap in Source / horizontal move)
        let mut curr_iy = neg_inf;
        
        // Boundary condition at j=0: Gap in Source (i monomers, 0 target)
        // This is valid for Ix state (vertical move) or M state
        let mut left_m = mode.gap_open + (i as i32) * mode.gap_extend;
        
        // Ix[i][0] should be consistent with M[i][0] at the boundary (pure gap)
        // We set it to the gap score so it can be extended in the next row.
        let curr_ix_0 = left_m;
        
        // Update prev arrays in-place as we scan? No, need old values.
        // Since we only need prev[j] and prev[j-1], we can iterate.
        // BUT simpler to allocate current row arrays to match `prev`.
        // Optimizing allocs: double buffering would be better, but let's just use a single curr buffer
        // or update prev carefully.
        // Let's use `curr_m` and `curr_ix` vectors for clarity, then swap.
        
        let mut curr_m_vec = vec![neg_inf; n + 1];
        let mut curr_ix_vec = vec![neg_inf; n + 1];
        
        curr_m_vec[0] = left_m;
        curr_ix_vec[0] = curr_ix_0;

        for j in 1..=n {
            // 1. Update Ix[i][j] (Vertical, gap in Target)
            // Depends on M[i-1][j] and Ix[i-1][j]
            let m_up = prev_m[j];
            let ix_up = prev_ix[j];
            let score_ix = (m_up + mode.gap_open).max(ix_up + mode.gap_extend);
            curr_ix_vec[j] = score_ix;

            // 2. Update Iy[i][j] (Horizontal, gap in Source)
            // Depends on M[i][j-1] and Iy[i][j-1]
            // left_m holds M[i][j-1] from previous col iteration
            // curr_iy holds Iy[i][j-1] from previous col iteration
            let m_left = left_m;
            let iy_left = curr_iy;
            let score_iy = (m_left + mode.gap_open).max(iy_left + mode.gap_extend);
            curr_iy = score_iy; // Update for next column

            // 3. Update M[i][j] (Match/Mismatch)
            // Depends on M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]
            // We need Iy[i-1][j-1]... wait.
            // Gotoh's recurrence:
            // M[i][j] = max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]) + score
            // We explicitly did NOT store Iy from the previous row.
            // However, Iy is "Insertion in Source".
            // Standard Gotoh definition:
            // M state usually collects from M, Ix, Iy.
            // But wait, if we defined Ix/Iy as "ending in gap", then:
            //   M[i][j] comes from (i-1, j-1). The state at (i-1, j-1) could be M, Ix, or Iy.
            //   Yes. So we DO need to know max(M, Ix, Iy) at (i-1, j-1).
            //   Let's define `prev_best[j]` ? 
            //   Or just re-calculate max at (i-1, j-1).
            //   Wait, standard implementation often simplifies:
            //   M[i][j] = Best[i-1][j-1] + score.
            //   Where Best[k][l] = max(M[k][l], Ix[k][l], Iy[k][l]).
            
            // Let's adjust storage: `prev_best` is enough if we don't distinguish source state for transition.
            // Actually, Gotoh usually implies:
            // Ix[i][j] = max(M[i-1][j] + open, Ix[i-1][j] + extend)
            // Iy[i][j] = max(M[i][j-1] + open, Iy[i][j-1] + extend)
            // M[i][j]  = max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]) + sub_score
            // So yes, we need the max value from (i-1, j-1).
            
            // Recover H[i-1][j-1] from `prev_best[j-1]`.
            let h_diag = prev_m[j-1]; // re-purposing prev_m as prev_best
            let score_match = h_diag + omega_monomer(&src[i-1], &tgt[j-1], mode);
            
            // Ix (Vertical) computation using `prev_best[j]`
            //   Ix[i][j] = max(Best[i-1][j] + open + extend, Ix[i-1][j] + extend)
            let h_up = prev_m[j];
            let ix_up = prev_ix[j];
            let score_ix = (h_up + mode.gap_open + mode.gap_extend).max(ix_up + mode.gap_extend);
            
            // Iy (Horizontal) computation
            //   Iy[i][j] = max(Best[i][j-1] + open + extend, Iy[i][j-1] + extend)
            //   Best[i][j-1] is `left_m` (curr_best of prev col)
            //   Iy[i][j-1] is `iy_left` (curr_iy of prev col) -- wait, iy_left is updated at step 2.
            //   Actually: Iy[i][j] = max(left_best + open + extend, iy_at_current + extend? No)
            //   We track `curr_iy` which represents Iy[i][j].
            //   curr_iy = max(left_best + open + extend, curr_iy_prev_col + extend)
            let score_iy = (left_m + mode.gap_open + mode.gap_extend).max(curr_iy + mode.gap_extend);
            
            // Final Best
            let best = score_match.max(score_ix).max(score_iy);
            
            // Store for next
            curr_m_vec[j] = best;
            curr_ix_vec[j] = score_ix;
            
            // Update running vars
            left_m = best;
            curr_iy = score_iy;
        }
        
        prev_m = curr_m_vec;
        prev_ix = curr_ix_vec;
    }
    prev_m[n]
}

pub struct HorCaches {
    pub match_at: HashMap<(usize, usize), MiniAlign>,
} // (a,b) -> internal alignment (optional)
impl HorCaches {
    pub fn new() -> Self {
        Self {
            match_at: HashMap::new(),
        }
    }
}

// Precise HOR tagging: for each monomer index, assign owning HOR id if within [start,end]
pub fn tag_monomers_by_hor(mons: &[MonRow], hors: &[HorRow]) -> Result<HorTagging, HorScanError> {
    // mons are assumed sorted by start; we still sort a local copy to be safe
    let mut sorted_mons = mons.to_vec();
    sorted_mons.sort_by_key(|r| r.start);
    let mut sorted_hors = hors.to_vec();
    sorted_hors.sort_by_key(|h| h.start);
    let mut units = Vec::new();
    let mut mon_to_hor = vec![None; sorted_mons.len() + 1]; // +1 for sentinel 0

    let mut k = 0usize; // index over monomers
    for (hid, h) in sorted_hors.iter().enumerate() {
        // find first mon fully inside HOR: mon.start >= h.start
        while k < sorted_mons.len() && sorted_mons[k].start < h.start {
            k += 1;
        }
        let start_k = k;
        // include while mon.end <= h.end (fully contained)
        while k < sorted_mons.len() && sorted_mons[k].end <= h.end {
            k += 1;
        }
        let end_k = k.saturating_sub(1);
        if start_k >= sorted_mons.len() || end_k < start_k {
            continue;
        }
        // store unit with same index
        let unit = HorUnit {
            label: h.hor_label.clone(),
            start_mi: start_k,
            end_mi: end_k,
        };
        for mi in unit.start_mi..=unit.end_mi {
            mon_to_hor[mi] = Some(hid);
        }
        units.push(unit);
    }
    Ok(HorTagging { mon_to_hor, units })
}

pub struct Incoming {
    pub match_map: HashMap<(usize, usize), Vec<(usize, usize, i32, usize, usize)>>, // (i,j)->[(si,sj,score,ai,bi)] sorted by (ai,bi)
    pub del_by_i: HashMap<usize, Vec<(usize, usize, i32)>>, // i -> [(si, hor_id, score)] sorted by hor_id
    pub ins_by_j: HashMap<usize, Vec<(usize, usize, i32)>>, // j -> [(sj, hor_id, score)] sorted by hor_id
}

pub fn build_incoming(
    src_units: &[HorUnit],
    tgt_units: &[HorUnit],
    pairs: &[(usize, usize, i32)],
    mode: &ModeHor,
) -> Incoming {
    let mut match_map: HashMap<(usize, usize), Vec<(usize, usize, i32, usize, usize)>> =
        HashMap::new();
    let mut del_by_i: HashMap<usize, Vec<(usize, usize, i32)>> = HashMap::new();
    let mut ins_by_j: HashMap<usize, Vec<(usize, usize, i32)>> = HashMap::new();

    for (ai, bi, score) in pairs {
        let a = &src_units[*ai];
        let b = &tgt_units[*bi];
        let i = a.end_mi + 1;
        let j = b.end_mi + 1;
        match_map
            .entry((i, j))
            .or_default()
            .push((a.start_mi, b.start_mi, *score, *ai, *bi));
    }
    for (hid, a) in src_units.iter().enumerate() {
        let i = a.end_mi + 1;
        let s = hor_gap(a.len(), mode);
        del_by_i.entry(i).or_default().push((a.start_mi, hid, s));
    }
    for (hid, b) in tgt_units.iter().enumerate() {
        let j = b.end_mi + 1;
        let s = hor_gap(b.len(), mode);
        ins_by_j.entry(j).or_default().push((b.start_mi, hid, s));
    }
    // Deterministic ordering for stable tie-break during backtrace
    for v in match_map.values_mut() {
        v.sort_by_key(|x| (x.3, x.4));
    }
    for v in del_by_i.values_mut() {
        v.sort_by_key(|x| x.1);
    }
    for v in ins_by_j.values_mut() {
        v.sort_by_key(|x| x.1);
    }
    let match_pairs: usize = match_map.values().map(|v| v.len()).sum();
    let del_count: usize = del_by_i.values().map(|v| v.len()).sum();
    let ins_count: usize = ins_by_j.values().map(|v| v.len()).sum();
    log::info!(
        "HOR transitions: matches={}, dels={}, ins={}",
        match_pairs,
        del_count,
        ins_count
    );
    Incoming {
        match_map,
        del_by_i,
        ins_by_j,
    }
}

pub fn nw_hierarchical(
    src: &[MonRow],
    tgt: &[MonRow],
    src_tag: &HorTagging,
    tgt_tag: &HorTagging,
    cache: &mut HorCaches,
    mode: &ModeHor,
    pairs: &[(usize, usize, i32)], // (a,b,HORScore)
    debug: bool,
) -> AlignOutputs {
    let m = src.len() - 1; // with sentinel
    let n = tgt.len() - 1;
    let cells_total = (m + 1) * (n + 1);
    // 3 x i32 (M/Ix/Iy) + 1 x u8 (Trace)
    let bytes_est = (cells_total as u64)
        * ((3 * std::mem::size_of::<i32>() + std::mem::size_of::<u8>()) as u64);
    let mb_est = (bytes_est as f64) / (1024.0 * 1024.0);
    log::info!(
        "DP matrices (affine): (m+1)x(n+1) = {}x{} (~{} cells), est RAM(M/Ix/Iy+Trace) ~{:.1} MiB",
        m + 1,
        n + 1,
        cells_total,
        mb_est
    );

    let cells = (m + 1) * (n + 1);
    let neg_inf = i32::MIN / 4;
    let go = mode.mon.gap_open;
    let ge = mode.mon.gap_extend;

    // Affine DP matrices
    let mut mat_m = vec![neg_inf; cells];
    let mut mat_ix = vec![neg_inf; cells];
    let mut mat_iy = vec![neg_inf; cells];
    let mut trace = vec![TraceCode::Empty as u8; cells]; // only for M state

    // init affine NW borders (global): gap(k) = go + k*ge
    mat_m[idx(0, 0, n)] = 0;
    // first column (j=0): prefixes of src vs empty tgt
    for i in 1..=m {
        let gap = go + (i as i32) * ge;
        let k = idx(i, 0, n);
        mat_m[k] = gap;
        mat_ix[k] = gap; // allow treating boundary as already in vertical gap
        // mat_iy[k] stays neg_inf
    }
    // first row (i=0): prefixes of tgt vs empty src
    for j in 1..=n {
        let gap = go + (j as i32) * ge;
        let k = idx(0, j, n);
        mat_m[k] = gap;
        mat_iy[k] = gap; // allow boundary horizontal gap
        // mat_ix[k] stays neg_inf
    }

    // incoming transitions (HOR-first priority later in tie-break)
    let incoming = build_incoming(&src_tag.units, &tgt_tag.units, pairs, mode);

    // We avoid storing per-cell debug steps to reduce memory; we'll emit during backtrace only.
    let mut dbg: Vec<DebugStep> = Vec::new();

    let stride = (m / 10).max(1); // ~10% progress granularity
    for i in 1..=m {
        if i % stride == 0 {
            log::info!("DP progress: row {}/{} (~{}%)", i, m, (i * 100) / m.max(1));
        }
        for j in 1..=n {
            let k = idx(i, j, n);
            let k_up = idx(i - 1, j, n);
            let k_left = idx(i, j - 1, n);
            let k_diag = idx(i - 1, j - 1, n);

            // 1) update Ix (vertical gap: source deletion / target insertion)
            let open_ix = mat_m[k_up].saturating_add(go).saturating_add(ge);
            let ext_ix = mat_ix[k_up].saturating_add(ge);
            mat_ix[k] = open_ix.max(ext_ix);

            // 2) update Iy (horizontal gap: source insertion / target deletion)
            let open_iy = mat_m[k_left].saturating_add(go).saturating_add(ge);
            let ext_iy = mat_iy[k_left].saturating_add(ge);
            mat_iy[k] = open_iy.max(ext_iy);

            // 3) micro alignment into M (match/mismatch)
            let score_base = omega_monomer(&src[i].mon, &tgt[j].mon, &mode.mon);
            let mut best_prev = mat_m[k_diag];
            let mut t_code = TraceCode::Diag; // from M
            if mat_ix[k_diag] > best_prev {
                best_prev = mat_ix[k_diag];
                t_code = TraceCode::FromIx;
            }
            if mat_iy[k_diag] > best_prev {
                best_prev = mat_iy[k_diag];
                t_code = TraceCode::FromIy;
            }
            let mut m_val = best_prev.saturating_add(score_base);

            // 4) HOR macro-jumps (MATCH / DEL / INS), all into M state with priority over micro
            // HOR MATCH: (si,sj,sc,ai,bi) with jump from (si,sj)->(i,j)
            if let Some(v) = incoming.match_map.get(&(i, j)) {
                for (si, sj, sc, _ai, _bi) in v {
                    let base = mat_m[idx(*si, *sj, n)];
                    let cand = base.saturating_add(*sc);
                    if cand > m_val || (cand == m_val && t_code != TraceCode::HorJump) {
                        m_val = cand;
                        t_code = TraceCode::HorJump;
                    }
                }
            }
            // HOR DEL: delete a HOR block in source (relative to target)
            if let Some(v) = incoming.del_by_i.get(&i) {
                for (si, _hid, sc) in v {
                    let base = mat_m[idx(*si, j, n)];
                    let cand = base.saturating_add(*sc);
                    if cand > m_val || (cand == m_val && t_code != TraceCode::HorJump) {
                        m_val = cand;
                        t_code = TraceCode::HorJump;
                    }
                }
            }
            // HOR INS: insert a HOR block in target (relative to source)
            if let Some(v) = incoming.ins_by_j.get(&j) {
                for (sj, _hid, sc) in v {
                    let base = mat_m[idx(i, *sj, n)];
                    let cand = base.saturating_add(*sc);
                    if cand > m_val || (cand == m_val && t_code != TraceCode::HorJump) {
                        m_val = cand;
                        t_code = TraceCode::HorJump;
                    }
                }
            }

            mat_m[k] = m_val;
            trace[k] = t_code as u8;
        }
    }

    // backtrace with macro expansion using precise HOR boundaries from tagging
    // State-aware traceback using M/Ix/Iy + Trace (HOR_JUMP + micro states)
    let mut rows = Vec::with_capacity(m + n);
    let mut hor_events: Vec<HorEvent> = Vec::new();
    let mut i = m;
    let mut j = n;
    // determine start state by best of M/Ix/Iy at (m,n)
    let end_k = idx(m, n, n);
    let mut best_end = mat_m[end_k];
    let mut state = 0u8; // 0: M, 1: Ix, 2: Iy
    if mat_ix[end_k] > best_end {
        best_end = mat_ix[end_k];
        state = 1;
    }
    if mat_iy[end_k] > best_end {
        // best_end = mat_iy[end_k]; // Unused
        state = 2;
    }
    let mut bt_steps: usize = 0;
    let bt_target_est = m + n; // rough upper bound (HOR expansions add rows in batches)
    let bt_stride = (bt_target_est / 20).max(10_000); // at most ~20 progress logs, min stride 10k
    while i > 0 || j > 0 {
        bt_steps += 1;
        if bt_steps % bt_stride == 0 {
            let pct = ((bt_steps.min(bt_target_est) * 100) / bt_target_est.max(1)).min(100);
            log::info!("Backtrace progress: ~{}% (steps={} i={} j={})", pct, bt_steps, i, j);
        }
        // State machine based traceback
        let k = idx(i, j, n);
        match state {
            0 => { // M state
                let code = match trace[k] {
                    x if x == TraceCode::HorJump as u8 => TraceCode::HorJump,
                    x if x == TraceCode::Diag as u8 => TraceCode::Diag,
                    x if x == TraceCode::FromIx as u8 => TraceCode::FromIx,
                    x if x == TraceCode::FromIy as u8 => TraceCode::FromIy,
                    _ => TraceCode::Empty,
                };
                match code {
                    TraceCode::HorJump => {
                        // Try HOR MATCH first
                        let mut took = false;
                        if let Some(v) = incoming.match_map.get(&(i, j)) {
                            for (si, sj, sc, ai, bi) in v {
                                let base = mat_m[idx(*si, *sj, n)];
                                if base.saturating_add(*sc) == mat_m[k] {
                                    if debug { dbg.push(DebugStep { i, j, step: StepCode::HorMatch }); }
                                    let au = &src_tag.units[*ai];
                                    let bu = &tgt_tag.units[*bi];
                                    let mini = if let Some(m) = cache.match_at.get(&(*ai, *bi)) { m.clone() } else {
                                        let src_slice: Vec<String> = (au.start_mi..=au.end_mi).map(|kk| src[kk + 1].mon.clone()).collect();
                                        let tgt_slice: Vec<String> = (bu.start_mi..=bu.end_mi).map(|kk| tgt[kk + 1].mon.clone()).collect();
                                        let mm = nw_global_mini(&src_slice, &tgt_slice, &mode.mon);
                                        cache.match_at.insert((*ai, *bi), mm.clone());
                                        mm
                                    };
                                    let mut si2 = au.start_mi + 1; // adjust for sentinel
                                    let mut tj2 = bu.start_mi + 1;
                                    let mut mth = 0usize; let mut mis = 0usize; let mut ins = 0usize; let mut del = 0usize;
                                    let mut diffs: Vec<i32> = Vec::new();
                                    let mut temp_rows: Vec<MonAlignmentRow> = Vec::new();
                                    for (a, b) in &mini.path {
                                        match (a, b) {
                                            (Some(_), Some(_)) => {
                                                let kind = if src[si2].mon == tgt[tj2].mon { mth += 1; "MTH" } else { mis += 1; "MIS" };
                                                let rel_src = (si2 - 1) as i32 - au.start_mi as i32;
                                                let rel_tgt = (tj2 - 1) as i32 - bu.start_mi as i32;
                                                diffs.push(rel_src - rel_tgt);
                                                push_mon_row(&mut temp_rows, &src[si2], &tgt[tj2], kind, si2, tj2);
                                                si2 += 1; tj2 += 1;
                                            }
                                            (Some(_), None) => { ins += 1; push_mon_row(&mut temp_rows, &src[si2], &tgt[0], "INS", si2, 0); si2 += 1; }
                                            (None, Some(_)) => { del += 1; push_mon_row(&mut temp_rows, &src[0], &tgt[tj2], "DEL", 0, tj2); tj2 += 1; }
                                            _ => {}
                                        }
                                    }
                                    temp_rows.reverse();
                                    rows.extend(temp_rows);
                                    let (mean_shift, shift_std, shift_consistent) = if diffs.is_empty() { (None, None, false) } else {
                                        let n_d = diffs.len() as f64; let sum: f64 = diffs.iter().map(|d| *d as f64).sum();
                                        let mean = sum / n_d; let var: f64 = diffs.iter().map(|d| { let x = *d as f64 - mean; x*x }).sum::<f64>() / n_d; let std = var.sqrt();
                                        let consistent = diffs.iter().all(|d| *d == diffs[0]) && diffs[0] != 0;
                                        (Some(mean), Some(std), consistent)
                                    };
                                    hor_events.push(HorEvent {
                                        event_type: "MATCH".into(),
                                        src_label: au.label.clone(),
                                        src_hor_id: Some(*ai),
                                        src_start_mi: au.start_mi,
                                        src_end_mi: au.end_mi,
                                        src_start_bp: src[au.start_mi + 1].start,
                                        src_end_bp: src[au.end_mi + 1].end,
                                        tgt_label: bu.label.clone(),
                                        tgt_hor_id: Some(*bi),
                                        tgt_start_mi: bu.start_mi,
                                        tgt_end_mi: bu.end_mi,
                                        tgt_start_bp: tgt[bu.start_mi + 1].start,
                                        tgt_end_bp: tgt[bu.end_mi + 1].end,
                                        score: mini.score + mode.hor_match_score,
                                        mth, mis, ins, del,
                                        mean_shift, shift_std, shift_consistent,
                                    });
                                    // Jump to the state BEFORE the HOR.
                                    // au.start_mi is 0-based index in src[1..].
                                    // DP index of start of HOR is au.start_mi + 1.
                                    // The state before that is DP index au.start_mi.
                                    i = au.start_mi;
                                    j = bu.start_mi;
                                    took = true;
                                    break;
                                }
                            }
                        }
                        if took { continue; }
                        // HOR DEL
                        if i > 0 {
                            if let Some(v) = incoming.del_by_i.get(&i) {
                                for (si, hid, sc) in v {
                                    let base = mat_m[idx(*si, j, n)];
                                    if base.saturating_add(*sc) == mat_m[k] {
                                        if debug { dbg.push(DebugStep { i, j, step: StepCode::HorDel }); }
                                        let unit = &src_tag.units[*hid];
                                        let mut temp_rows: Vec<MonAlignmentRow> = Vec::new();
                                        for kk in unit.start_mi..=unit.end_mi {
                                            let kk2 = kk + 1;
                                            push_mon_row(&mut temp_rows, &src[kk2], &tgt[0], "INS", kk2, 0);
                                        }
                                        temp_rows.reverse(); rows.extend(temp_rows);
                                        hor_events.push(HorEvent {
                                            event_type: "INS".into(),
                                            src_label: unit.label.clone(),
                                            src_hor_id: Some(*hid),
                                            src_start_mi: unit.start_mi,
                                            src_end_mi: unit.end_mi,
                                            src_start_bp: src[unit.start_mi + 1].start,
                                            src_end_bp: src[unit.end_mi + 1].end,
                                            tgt_label: "-".into(),
                                            tgt_hor_id: None,
                                            tgt_start_mi: 0,
                                            tgt_end_mi: 0,
                                            tgt_start_bp: 0,
                                            tgt_end_bp: 0,
                                            score: hor_gap(unit.len(), mode),
                                            mth: 0, mis: 0, ins: unit.len(), del: 0,
                                            mean_shift: None, shift_std: None, shift_consistent: false,
                                        });
                                        i = unit.start_mi;
                                        took = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if took { continue; }
                        // HOR INS
                        if j > 0 {
                            if let Some(v) = incoming.ins_by_j.get(&j) {
                                for (sj, hid, sc) in v {
                                    let base = mat_m[idx(i, *sj, n)];
                                    if base.saturating_add(*sc) == mat_m[k] {
                                        if debug { dbg.push(DebugStep { i, j, step: StepCode::HorIns }); }
                                        let unit = &tgt_tag.units[*hid];
                                        let mut temp_rows: Vec<MonAlignmentRow> = Vec::new();
                                        for kk in unit.start_mi..=unit.end_mi {
                                            let kk2 = kk + 1;
                                            push_mon_row(&mut temp_rows, &src[0], &tgt[kk2], "DEL", 0, kk2);
                                        }
                                        temp_rows.reverse(); rows.extend(temp_rows);
                                        hor_events.push(HorEvent {
                                            event_type: "DEL".into(),
                                            src_label: "-".into(),
                                            src_hor_id: None,
                                            src_start_mi: 0,
                                            src_end_mi: 0,
                                            src_start_bp: 0,
                                            src_end_bp: 0,
                                            tgt_label: unit.label.clone(),
                                            tgt_hor_id: Some(*hid),
                                            tgt_start_mi: unit.start_mi,
                                            tgt_end_mi: unit.end_mi,
                                            tgt_start_bp: tgt[unit.start_mi + 1].start,
                                            tgt_end_bp: tgt[unit.end_mi + 1].end,
                                            score: hor_gap(unit.len(), mode),
                                            mth: 0, mis: 0, ins: 0, del: unit.len(),
                                            mean_shift: None, shift_std: None, shift_consistent: false,
                                        });
                                        j = unit.start_mi;
                                        took = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if took { continue; }
                        panic!("Backtrace HOR_JUMP unresolved at ({},{})", i, j);
                    }
                    TraceCode::Diag | TraceCode::FromIx | TraceCode::FromIy => {
                        // micro alignment step: (i-1,j-1) -> (i,j)
                        if i == 0 || j == 0 {
                            panic!("Micro step at boundary ({},{})", i, j);
                        }
                        let kind = if src[i].mon == tgt[j].mon { "MTH" } else { "MIS" };
                        if debug {
                            dbg.push(DebugStep { i, j, step: StepCode::Diag });
                        }
                        push_mon_row(&mut rows, &src[i], &tgt[j], kind, i, j);
                        i -= 1;
                        j -= 1;
                        state = match code {
                            TraceCode::Diag => 0,
                            TraceCode::FromIx => 1,
                            TraceCode::FromIy => 2,
                            _ => 0,
                        };
                        continue;
                    }
                    TraceCode::Empty => {
                        // should only happen at boundaries; fall through to gap logic below
                    }
                }
            }
            1 => { // Ix state: vertical gap
                if i == 0 {
                    panic!("Ix state at top row j={} with i=0", j);
                }
                let k_up = idx(i - 1, j, n);
                let curr = mat_ix[k];
                let open_ix = mat_m[k_up].saturating_add(go).saturating_add(ge);
                let ext_ix = mat_ix[k_up].saturating_add(ge);
                if curr == ext_ix {
                    if debug { dbg.push(DebugStep { i, j, step: StepCode::Up }); }
                    push_mon_row(&mut rows, &src[i], &tgt[0], "INS", i, 0);
                    i -= 1;
                    state = 1;
                    continue;
                } else if curr == open_ix {
                    if debug { dbg.push(DebugStep { i, j, step: StepCode::Up }); }
                    push_mon_row(&mut rows, &src[i], &tgt[0], "INS", i, 0);
                    i -= 1;
                    state = 0;
                    continue;
                } else {
                    panic!("Ix traceback mismatch at ({},{})", i, j);
                }
            }
            2 => { // Iy state: horizontal gap
                if j == 0 {
                    panic!("Iy state at left col i={} with j=0", i);
                }
                let k_left = idx(i, j - 1, n);
                let curr = mat_iy[k];
                let open_iy = mat_m[k_left].saturating_add(go).saturating_add(ge);
                let ext_iy = mat_iy[k_left].saturating_add(ge);
                if curr == ext_iy {
                    if debug { dbg.push(DebugStep { i, j, step: StepCode::Left }); }
                    push_mon_row(&mut rows, &src[0], &tgt[j], "DEL", 0, j);
                    j -= 1;
                    state = 2;
                    continue;
                } else if curr == open_iy {
                    if debug { dbg.push(DebugStep { i, j, step: StepCode::Left }); }
                    push_mon_row(&mut rows, &src[0], &tgt[j], "DEL", 0, j);
                    j -= 1;
                    state = 0;
                    continue;
                } else {
                    panic!("Iy traceback mismatch at ({},{})", i, j);
                }
            }
            _ => unreachable!(),
        }
        // If we reach here with i==0 or j==0, fall back to consuming remaining gaps in chosen state
        if i > 0 && j == 0 {
            // must be vertical gaps
            push_mon_row(&mut rows, &src[i], &tgt[0], "INS", i, 0);
            i -= 1;
            state = 1;
            continue;
        }
        if j > 0 && i == 0 {
            push_mon_row(&mut rows, &src[0], &tgt[j], "DEL", 0, j);
            j -= 1;
            state = 2;
            continue;
        }
        if i == 0 && j == 0 { break; }
        panic!("Backtrace fell through at ({},{}) state={}", i, j, state);
    }
    rows.reverse();
    hor_events.reverse(); // events were collected from end to start
    log::info!("Backtrace done: steps={} rows={} hor.events={}", bt_steps, rows.len(), hor_events.len());

    // ============== Consistency validation (coverage & HOR event accounting) ==============
    {
        // Coverage arrays: DP indices 1..=m / 1..=n correspond to real monomers
        let mut src_cov = vec![0usize; m + 1];
        let mut tgt_cov = vec![0usize; n + 1];
        for r in &rows {
            if r.src_mi > 0 && r.src_mon != "-" && r.src_mi <= m { src_cov[r.src_mi] += 1; }
            if r.tgt_mi > 0 && r.tgt_mon != "-" && r.tgt_mi <= n { tgt_cov[r.tgt_mi] += 1; }
        }
        let mut src_ok = true; let mut tgt_ok = true;
        for i2 in 1..=m { if src_cov[i2] != 1 { src_ok = false; log::warn!("Source monomer DP index {} seen {} times (expect 1)", i2, src_cov[i2]); } }
        for j2 in 1..=n { if tgt_cov[j2] != 1 { tgt_ok = false; log::warn!("Target monomer DP index {} seen {} times (expect 1)", j2, tgt_cov[j2]); } }
        if src_ok && tgt_ok { log::info!("Monomer coverage validation passed: all non-gap monomers appear exactly once."); }
        else { log::warn!("Monomer coverage validation FAILED (see warnings above). Potential coordinate confusion."); }

        // HOR INS / DEL event internal consistency
        for e in &hor_events {
            match e.event_type.as_str() {
                "INS" => {
                    let len = e.src_end_mi.saturating_sub(e.src_start_mi) + 1;
                    if e.ins != len { log::warn!("HOR INS event {} length {} but ins count {}", e.src_label, len, e.ins); }
                }
                "DEL" => {
                    let len = e.tgt_end_mi.saturating_sub(e.tgt_start_mi) + 1;
                    if e.del != len { log::warn!("HOR DEL event {} length {} but del count {}", e.tgt_label, len, e.del); }
                }
                _ => {}
            }
        }
    }
    // ================================================================================

    AlignOutputs {
        rows,
        debug: dbg,
        hor_events,
    }
}

fn push_mon_row(out: &mut Vec<MonAlignmentRow>, s: &MonRow, t: &MonRow, kind: &str, src_mi: usize, tgt_mi: usize) {
    let (src_row, tgt_row, smi, tmi) = if kind == "DEL" { // source deletion relative to target present (source gap)
        (
            MonRow { sample: s.sample.clone(), start: s.start, end: s.start, mon: "-".into() },
            t.clone(),
            0usize,
            tgt_mi,
        )
    } else if kind == "INS" { // source insertion relative to target (target gap)
        (
            s.clone(),
            MonRow { sample: t.sample.clone(), start: t.start, end: t.start, mon: "-".into() },
            src_mi,
            0usize,
        )
    } else { // MTH or MIS
        (s.clone(), t.clone(), src_mi, tgt_mi)
    };
    out.push(MonAlignmentRow {
        src_sample: src_row.sample,
        src_start: src_row.start,
        src_end: src_row.end,
        src_mon: src_row.mon,
        tgt_sample: tgt_row.sample,
        tgt_start: tgt_row.start,
        tgt_end: tgt_row.end,
        tgt_mon: tgt_row.mon,
        align_type: kind.into(),
        src_mi: smi,
        tgt_mi: tmi,
    });
}
