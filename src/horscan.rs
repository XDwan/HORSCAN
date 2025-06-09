use crate::io::{read_bed_file, HorScanError, Mode, MonRow};
use crate::optimize::optimize_blast_file;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn horscan_main(
    source_bed_path: String,
    target_bed_path: String,
    output_prefix: String,
    mode: Mode,
) -> Result<(), HorScanError> {
    // 加载hor文件
    let path = output_prefix.clone() + ".alignment";
    // match fs::metadata(path.clone()) {
    //     Ok(_) => {
    //         // 如果文件存在，则终止程序
    //         println!("File {} exists, terminating program.", path.clone());
    //         process::exit(0); // 退出程序，0 表示正常退出
    //     }
    //     Err(_) => {}
    // }

    let mut source_bed = read_bed_file(source_bed_path.as_str())?;
    let mut target_bed = read_bed_file(target_bed_path.as_str())?;
    source_bed.sort_by_key(|row| row.start);
    target_bed.sort_by_key(|row| row.start);
    let source_sample = source_bed[0].sample.clone();
    let target_sample = target_bed[0].sample.clone();
    print!(
        "source sample: {} Monomer Length: {}\n",
        source_sample,
        source_bed.len()
    );
    print!(
        "target sample: {} Monomer Length: {}\n",
        target_sample,
        target_bed.len()
    );
    source_bed.insert(
        0,
        MonRow {
            sample: source_sample.clone(),
            start: 0,
            end: 0,
            mon: "-".to_string(),
        },
    );
    target_bed.insert(
        0,
        MonRow {
            sample: target_sample.clone(),
            start: 0,
            end: 0,
            mon: "-".to_string(),
        },
    );
    // 调用HORSCAN的动态规划算法

    println!("memory update linear score");
    let score = horscan_global_linear_score_alignment_low(
        source_bed.clone(),
        target_bed.clone(),
        mode.clone(),
    )?;
    let alignment_path = find_global_alignment_all_path_low(&score)?;
    // 保存文件
    save_all_path(
        path.clone(),
        source_bed.clone(),
        target_bed.clone(),
        &alignment_path,
    )?;

    // optimize alignment path
    optimize_blast_file(path.clone())?;

    // 保存 score matrix
    // save_score_matrix(&score, path.clone())?;
    Ok(())
}

fn horscan_global_linear_score_alignment_low(
    source_bed: Vec<MonRow>,
    target_bed: Vec<MonRow>,
    mode: Mode,
) -> Result<Vec<Vec<i32>>, HorScanError> {
    let source_len = source_bed.len();
    let target_len = target_bed.len();
    let mut score_matrix = vec![vec![0; target_len]; source_len];
    for i in 0..source_len {
        score_matrix[i][0] = mode.gap_score * i as i32;
    }
    for j in 0..target_len {
        score_matrix[0][j] = mode.gap_score * j as i32;
    }

    // 将连续的q或者t gap 的gap进行压缩
    let mut last_path = 0;

    for i in 1..source_len {
        for j in 1..target_len {
            let match_score = if source_bed[i].mon == target_bed[j].mon {
                score_matrix[i - 1][j - 1] + mode.match_score
            } else {
                score_matrix[i - 1][j - 1] + mode.mismatch_score
            };
            let mut q_gap_score = score_matrix[i - 1][j] + mode.gap_score;
            let mut t_gap_socre = score_matrix[i][j - 1] + mode.gap_score;
            // if last_path == 1 {
            //     q_gap_score = score_matrix[i - 1][j] + mode.gap_score / 2;
            // } else if last_path == 2 {
            //     t_gap_socre = score_matrix[i][j - 1] + mode.gap_score / 2;
            // }
            // 如果当前的gap是连续的，那么就添加half gap score
            let last_score = vec![q_gap_score, t_gap_socre, match_score, 0];

            let max_score = *last_score.iter().max().unwrap();
            // last_path = last_score.iter().position(|&x| x == max_score).unwrap();

            score_matrix[i][j] = max_score;
        }
    }
    Ok(score_matrix)
}

fn find_global_alignment_all_path_low(
    score_matrix: &Vec<Vec<i32>>,
) -> Result<Vec<(usize, usize, usize)>, HorScanError> {
    // find the path with the highest score
    let mut i = score_matrix.len() - 1;
    let mut j = score_matrix[0].len() - 1;
    let mut path = Vec::new();
    // 0: MTH/MIS 1: INS 2:DEL
    print!("source length: {} target length: {}\n", i, j);

    // init: i j  find max score in  socre_matrix[i][:] and score_matrix[:][j]
    let mut max_score = score_matrix[i][j]; // start from the bottom-right corner
    let mut max_i = i;
    let mut max_j = j;
    for k in 0..i {
        if score_matrix[k][j] > max_score {
            max_score = score_matrix[k][j];
            max_i = k;
        }
    }
    for k in 0..j {
        if score_matrix[i][k] > max_score {
            max_score = score_matrix[i][k];
            max_j = k;
        }
    }
    // add break path
    if max_i == i && max_j == j {
        // path.push((i, j, 0 as usize));
    } else if score_matrix[max_i][j] > score_matrix[i][max_j] {
        for k in (max_i + 1..=i).rev() {
            path.push((k, j, 1 as usize));
        }
        i = max_i;
    } else {
        for k in (max_j + 1..=j).rev() {
            path.push((i, k, 2 as usize));
        }
        j = max_j;
    }
    path.push((i, j, 0 as usize));
    while i > 0 && j > 0 {
        let last_path: Vec<(usize, usize)> = vec![(i - 1, j - 1), (i - 1, j), (i, j - 1)];
        let last_score: Vec<i32> = vec![
            score_matrix[i - 1][j - 1],
            score_matrix[i - 1][j],
            score_matrix[i][j - 1],
        ];
        let max_score = last_score.iter().max().unwrap();
        let max_index = last_score.iter().position(|&x| x == *max_score).unwrap();
        i = last_path[max_index].0;
        j = last_path[max_index].1;
        if i > 0 && j > 0 {
            path.push((i, j, max_index as usize));
        }
    }
    if i > 0 {
        for k in (0..=i).rev() {
            path.push((k, 0, 1));
        }
    } else if j > 0 {
        for k in (0..=j).rev() {
            path.push((0, k, 2));
        }
    }

    Ok(path)
}

fn save_all_path(
    file_path: String,
    source_bed: Vec<MonRow>,
    target_bed: Vec<MonRow>,
    alignment_path: &Vec<(usize, usize, usize)>,
) -> Result<(), HorScanError> {
    let mut file = File::create(file_path).expect("无法创建文件");
    // 0: MTH/MIS 1: INS 2:DEL
    let alignment_type = vec!["MTH", "INS", "DEL", "MIS"];
    for &path in alignment_path.iter().rev() {
        let mut source_row = source_bed[path.0].clone();
        let mut target_row = target_bed[path.1].clone();
        let mut alignment_type = alignment_type[path.2];
        if path.0 == 0 && path.1 == 0 {
            continue;
        }
        if path.2 == 2 {
            source_row = MonRow {
                sample: source_row.sample.clone(),
                start: source_row.start,
                end: source_row.start,
                mon: "-".to_string(),
            };
        } else if path.2 == 1 {
            target_row = MonRow {
                sample: target_row.sample.clone(),
                start: target_row.start,
                end: target_row.start,
                mon: "-".to_string(),
            };
        } else if path.2 == 0 {
            if source_row.mon == target_row.mon {
                alignment_type = "MTH";
            } else {
                alignment_type = "MIS";
            }
        } else {
            // 抛出错误 AlignmentError
            return Err(HorScanError::AlignmentError("未知的对齐类型".to_string()));
        }
        // 输出模式
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            source_row.sample,
            source_row.start,
            source_row.end,
            source_row.mon,
            target_row.sample,
            target_row.start,
            target_row.end,
            target_row.mon,
            alignment_type,
        )?;
        // 调试模式
        // writeln!(
        //     file,
        //     "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        //     source_row.sample,
        //     source_row.start,
        //     source_row.end,
        //     source_row.mon,
        //     target_row.sample,
        //     target_row.start,
        //     target_row.end,
        //     target_row.mon,
        //     alignment_type,
        //     path.0,
        //     path.1,
        //     path.2,
        // )?;
    }
    Ok(())
}

// 将 score 保存为numpy能够直接打开的格式
fn save_score_matrix(
    score_matrix: &Vec<Vec<i32>>,
    output_path: String,
) -> Result<(), HorScanError> {
    let file = File::create(output_path + ".mat")?;
    let mut writer = BufWriter::new(file);
    for row in score_matrix {
        for score in row {
            write!(writer, "{} ", score)?; // 使用空格分隔每个分数
        }
        writeln!(writer)?;
    }
    Ok(())
}
