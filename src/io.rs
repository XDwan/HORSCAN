use std::{fs::File, io::{BufRead, BufReader, Write, BufWriter}, path::Path};
use crate::{errors::HorScanError, types::*};

pub fn read_mon_bed<P: AsRef<Path>>(p: P) -> Result<Vec<MonRow>, HorScanError> {
    let f = File::open(&p).map_err(|e| HorScanError::Io(format!("{}: {e}", p.as_ref().display())))?;
    let rdr = BufReader::new(f);
    let mut out = Vec::new();
    for (ln, line) in rdr.lines().enumerate() {
        let line = line.map_err(|e| HorScanError::Io(format!("{}: {e}", p.as_ref().display())))?;
        if line.trim().is_empty() || line.starts_with('#') || line.starts_with("track") { continue; }
        let cols: Vec<_> = line.split_whitespace().collect();
        if cols.len() < 4 {
            return Err(HorScanError::ParseError(format!("{}:{} need 4 cols", p.as_ref().display(), ln+1)));
        }
        let row = MonRow {
            sample: cols[0].to_string(),
            start: cols[1].parse().map_err(|_| HorScanError::ParseError(format!("bad start @{}:{}", p.as_ref().display(), ln+1)))?,
            end: cols[2].parse().map_err(|_| HorScanError::ParseError(format!("bad end @{}:{}", p.as_ref().display(), ln+1)))?,
            mon: cols[3].to_string(),
        };
        out.push(row);
    }
    Ok(out)
}

pub fn read_hor_bed<P: AsRef<Path>>(p: P) -> Result<Vec<HorRow>, HorScanError> {
    let f = File::open(&p).map_err(|e| HorScanError::Io(format!("{}: {e}", p.as_ref().display())))?;
    let rdr = BufReader::new(f);
    let mut out = Vec::new();
    for (ln, line) in rdr.lines().enumerate() {
        let line = line.map_err(|e| HorScanError::Io(format!("{}: {e}", p.as_ref().display())))?;
        if line.trim().is_empty() || line.starts_with('#') || line.starts_with("track") { continue; }
        let cols: Vec<_> = line.split_whitespace().collect();
        if cols.len() < 4 {
            return Err(HorScanError::ParseError(format!("{}:{} need 4 cols", p.as_ref().display(), ln+1)));
        }
        let row = HorRow {
            sample: cols[0].to_string(),
            start: cols[1].parse().map_err(|_| HorScanError::ParseError(format!("bad start @{}:{}", p.as_ref().display(), ln+1)))?,
            end: cols[2].parse().map_err(|_| HorScanError::ParseError(format!("bad end @{}:{}", p.as_ref().display(), ln+1)))?,
            hor_label: cols[3].to_string(),
        };
        out.push(row);
    }
    Ok(out)
}

// Save only the primary alignment (logical order). Debug matrix output removed to reduce I/O.
pub fn save_alignment_outputs<P: AsRef<Path>>(prefix: P, out: &AlignOutputs) -> Result<(), HorScanError> {
    let align_path = {
        let p = prefix.as_ref();
        let mut out = p.to_path_buf();
        if let Some(name) = p.file_name() {
            let mut new_name = name.to_os_string();
            new_name.push(".alignment.tsv");
            out.set_file_name(new_name);
        } else {
            out.push(".alignment.tsv");
        }
        out
    };
    let f = File::create(&align_path)
        .map_err(|e| HorScanError::Io(format!("{}: {e}", align_path.display())))?;
    let mut w = BufWriter::new(f);
    for r in &out.rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.src_sample, r.src_start, r.src_end, r.src_mon,
            r.tgt_sample, r.tgt_start, r.tgt_end, r.tgt_mon,
            r.align_type, r.src_mi, r.tgt_mi
        ).map_err(|e| HorScanError::Io(format!("{}: {e}", align_path.display())))?;
    }
    Ok(())
}

pub fn save_hor_events<P: AsRef<Path>>(prefix: P, events: &[HorEvent]) -> Result<(), HorScanError> {
    let path = {
        let p = prefix.as_ref();
        let mut out = p.to_path_buf();
        if let Some(name) = p.file_name() {
            let mut new_name = name.to_os_string();
            new_name.push(".hor.tsv");
            out.set_file_name(new_name);
        } else {
            out.push(".hor.tsv");
        }
        out
    };
    let f = File::create(&path)
        .map_err(|e| HorScanError::Io(format!("{}: {e}", path.display())))?;
    let mut w = BufWriter::new(f);
    // header
    writeln!(w, "event\tsrc.label\tsrc.hor.id\tsrc.mi.start\tsrc.mi.end\tsrc.bp.start\tsrc.bp.end\ttgt.label\ttgt.hor.id\ttgt.mi.start\ttgt.mi.end\ttgt.bp.start\ttgt.bp.end\tscore\tmth\tmis\tins\tdel\tmean.shift\tshift.std\tshift.consistent")
        .map_err(|e| HorScanError::Io(format!("{}: {e}", path.display())))?;
    for e in events {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.event_type,
            e.src_label, e.src_hor_id.map(|v| v.to_string()).unwrap_or("NA".into()), e.src_start_mi, e.src_end_mi, e.src_start_bp, e.src_end_bp,
            e.tgt_label, e.tgt_hor_id.map(|v| v.to_string()).unwrap_or("NA".into()), e.tgt_start_mi, e.tgt_end_mi, e.tgt_start_bp, e.tgt_end_bp,
            e.score, e.mth, e.mis, e.ins, e.del,
            e.mean_shift.map(|v| v.to_string()).unwrap_or("NA".into()),
            e.shift_std.map(|v| v.to_string()).unwrap_or("NA".into()),
            if e.shift_consistent { "1" } else { "0" }
        ).map_err(|e| HorScanError::Io(format!("{}: {e}", path.display())))?;
    }
    Ok(())
}

// Removed: save_stats, save_shift_rows, save_unaligned, save_alignment_sorted to minimize outputs.
