use crate::{
    align::{nw_hierarchical, tag_monomers_by_hor, HorCaches},
    args::{validate_args, Args},
    errors::HorScanError,
    io::{read_hor_bed, read_mon_bed},
    types::*,
};
use anyhow::Result;

pub fn horscan_main(cfg: Args) -> Result<()> {
    validate_args(&cfg)?;

    log::info!("HORSCANv starting...\n source mon: {}\n target mon: {}\n source HOR: {}\n target HOR: {}\n output: {}",
        cfg.source, cfg.target, cfg.source_hor, cfg.target_hor, cfg.output_prefix);

    // Normalize / ensure output directory exists (support running from repo root OR subdir)
    let out_prefix = std::path::PathBuf::from(&cfg.output_prefix);
    if let Some(parent) = out_prefix.parent() {
        if !parent.exists() {
            if let Err(e) = std::fs::create_dir_all(parent) {
                return Err(anyhow::anyhow!(HorScanError::Io(format!("create {}: {e}", parent.display()))));
            }
        }
    }

    // Convert magnitudes to signed scores (Affine only): mismatch/gap/macro mismatch are penalties
    // Micro (monomer) set
    let match_score = cfg.mon_scores[0];
    let mismatch_score = -cfg.mon_scores[1].abs();
    let gap_open = -cfg.mon_scores[2].abs();
    let gap_extend = -cfg.mon_scores[3].abs();

    // Macro (HOR) set
    let hor_match_score = cfg.hor_scores[0];
    let hor_mismatch_score = -cfg.hor_scores[1].abs();
    let hor_gap_open = -cfg.hor_scores[2].abs();
    let hor_gap_extend = -cfg.hor_scores[3].abs();

    let mode = ModeHor {
        mon: ModeMon {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        },
        relaxed: cfg.relaxed,
        hor_match_score,
        hor_mismatch_score,
        hor_gap_open,
        hor_gap_extend,
    };

    // Read monomers
    let mut src = read_mon_bed(&cfg.source)?;
    let mut tgt = read_mon_bed(&cfg.target)?;
    src.sort_by_key(|r| r.start);
    tgt.sort_by_key(|r| r.start);
    let src_sentinel = MonRow {
        sample: src
            .first()
            .map(|x| x.sample.clone())
            .unwrap_or_else(|| "SRC".into()),
        start: 0,
        end: 0,
        mon: "-".into(),
    };
    let tgt_sentinel = MonRow {
        sample: tgt
            .first()
            .map(|x| x.sample.clone())
            .unwrap_or_else(|| "TGT".into()),
        start: 0,
        end: 0,
        mon: "-".into(),
    };
    src.insert(0, src_sentinel);
    tgt.insert(0, tgt_sentinel);

    // Read HORs and tag monomers by HOR
    let mut src_h = read_hor_bed(&cfg.source_hor)?;
    let mut tgt_h = read_hor_bed(&cfg.target_hor)?;
    src_h.sort_by_key(|r| r.start);
    tgt_h.sort_by_key(|r| r.start);

    // Tagging uses real mon rows without sentinel
    let src_tag = tag_monomers_by_hor(&src[1..], &src_h)?;
    let tgt_tag = tag_monomers_by_hor(&tgt[1..], &tgt_h)?;
    log::info!(
        "Loaded monomers: src={} tgt={}",
        src.len() - 1,
        tgt.len() - 1
    );
    log::info!("Loaded HORs: src={} tgt={}", src_h.len(), tgt_h.len());

    // Ω strict: same label only
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let band = cfg.hor_pair_band;
    for (ai, a) in src_tag.units.iter().enumerate() {
        for (bi, b) in tgt_tag.units.iter().enumerate() {
            if a.label != b.label {
                continue;
            }
            if band > 0 {
                // band heuristic: constrain absolute difference of end monomer index (coarse colinearity)
                let diff = if a.end_mi > b.end_mi {
                    a.end_mi - b.end_mi
                } else {
                    b.end_mi - a.end_mi
                };
                if diff > band {
                    continue;
                }
            }
            pairs.push((ai, bi));
        }
    }

    // Precompute HOR-MATCH scores and cache mini alignments
    let mut cache = HorCaches::new();
    let mut pairs_scored: Vec<(usize, usize, i32)> = Vec::with_capacity(pairs.len());
    for (k, (ai, bi)) in pairs.into_iter().enumerate() {
        let a = &src_tag.units[ai];
        let b = &tgt_tag.units[bi];
        let src_slice: Vec<String> = (a.start_mi..=a.end_mi)
            .map(|k| src[k + 1].mon.clone()) // +1 to skip sentinel
            .collect();
        let tgt_slice: Vec<String> = (b.start_mi..=b.end_mi)
            .map(|k| tgt[k + 1].mon.clone())
            .collect();
        // Micro sum + Macro match bonus
        let micro_score = crate::align::nw_global_mini_score(&src_slice, &tgt_slice, &mode.mon);
        let total = micro_score + mode.hor_match_score;
        pairs_scored.push((ai, bi, total));
        if (k + 1) % 100000 == 0 {
            log::info!(
                "Precomputed {} HOR pairs ({}%)...",
                k + 1,
                ((k + 1) * 100) / pairs_scored.capacity().max(1)
            );
        }
    }

    log::info!("Scoring params: micro(match={}, mismatch={}, open={}, extend={}), macro(match={}, mismatch={}, open={}, extend={})",
        mode.mon.match_score, mode.mon.mismatch_score, mode.mon.gap_open, mode.mon.gap_extend,
        mode.hor_match_score, mode.hor_mismatch_score, mode.hor_gap_open, mode.hor_gap_extend);
    let outputs = nw_hierarchical(
        &src,
        &tgt,
        &src_tag,
        &tgt_tag,
        &mut cache,
        &mode,
        &pairs_scored,
        cfg.debug,
    );

    // Sort HOR events by target then source coordinate for deterministic downstream analysis
    let mut events_sorted = outputs.hor_events.clone();
    events_sorted.sort_by_key(|e| (e.tgt_start_bp, e.src_start_bp));
    crate::io::save_hor_events(&out_prefix, &events_sorted)?;
    crate::io::save_alignment_outputs(&out_prefix, &outputs)?;

    log::info!(
        "Done. Alignment rows: {}, HOR events: {}",
        outputs.rows.len(),
        events_sorted.len()
    );
    Ok(())
}
