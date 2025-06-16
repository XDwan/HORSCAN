use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::io::{BufRead, BufReader};
use std::path::Path;
use thiserror::Error; // 引入 thiserror 的宏

#[derive(Clone, Debug)]
pub struct Mode {
    // pub match_method: i32,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_score: i32,
    // pub jump_max: usize,
}

#[derive(Clone, Debug)]
pub struct MonRow {
    pub sample: String,
    pub start: i32,
    pub end: i32,
    pub mon: String,
}

#[derive(Clone, Debug)]
pub struct MonAlignmentRow {
    pub source: String,
    pub source_start: i32,
    pub source_end: i32,
    pub source_mon: String,
    pub target: String,
    pub target_start: i32,
    pub target_end: i32,
    pub target_mon: String,
    pub align_type: String
}

// 修正错误枚举定义
#[derive(Debug, Error)]
pub enum HorScanError {
    // 文件操作错误（自动转换 std::io::Error）
    #[error("IO error: {0}")]
    Io(#[from] io::Error),

    // BED文件解析错误
    #[error("Parse error: {0}")]
    ParseError(String),

    // 动态规划逻辑错误
    #[error("Alignment error: {0}")]
    AlignmentError(String)
}

/// 辅助函数：将字符串解析为特定类型，处理转换错误
fn parse_field<T: std::str::FromStr>(s: &str, field_name: &str) -> Result<T, HorScanError> {
    s.parse()
        .map_err(|_| HorScanError::ParseError(format!("无法解析{}字段: '{}'", field_name, s)))
}

/// 读取BED文件并转换为MonRow结构体
pub fn read_bed_file<P: AsRef<Path>>(path: P) -> Result<Vec<MonRow>, HorScanError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        // 跳过注释行和空行
        if line.starts_with("track") || line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            return Err(HorScanError::ParseError(format!(
                "行 '{}' 需要至少4列",
                line
            )));
        }

        // 解析前三个字段为i32
        let sample = parts[0].to_string();
        let start = parse_field(parts[1], "start")?;
        let end = parse_field(parts[2], "end")?;

        // 第四个字段保持String类型
        let mon = parts[3].to_string();

        rows.push(MonRow {
            sample,
            start,
            end,
            mon,
        });
    }
    Ok(rows)
}

pub fn read_blast_file<P: AsRef<Path>>(path: P) -> Result<Vec<MonAlignmentRow>, HorScanError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        // 跳过注释行和空行
        if line.starts_with("track") || line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            return Err(HorScanError::ParseError(format!(
                "行 '{}' 需要至少9列",
                line
            )));
        }

        // 解析前三个字段为i32
        let source = parts[0].to_string();
        let source_start = parse_field(parts[1], "start")?;
        let source_end = parse_field(parts[2], "end")?;
        let source_mon = parts[3].to_string();
        let target = parts[4].to_string();
        let target_start = parse_field(parts[5], "start")?;
        let target_end = parse_field(parts[6], "end")?;
        let target_mon = parts[7].to_string();
        let align_type = parts[8].to_string();

        // 第四个字段保持String类型
        // let mon = parts[3].to_string();

        rows.push(MonAlignmentRow {
            source,
            source_start,
            source_end,
            source_mon,
            target,
            target_start,
            target_end,
            target_mon,
            align_type,
        });
    }
    Ok(rows)
}

pub fn save_blast_file(rows: &Vec<MonAlignmentRow>, output_path: String) -> Result<(), HorScanError> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    for row in rows {
        writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                 row.source, row.source_start, row.source_end, row.source_mon, 
                 row.target, row.target_start, row.target_end, row.target_mon, 
                 row.align_type)?;

    }
    Ok(())

}
// 测试read_bed函数
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_read_bed() -> Result<(), HorScanError> {
        let bed_path = "/home/wan/projects/HORSCAN/test_output/G100_monomer.bed";
        let result = read_bed_file(bed_path)?;
        // 逐行打印
        for row in result {
            println!("{:?}\n", row);
        }
        Ok(())
    }
}
