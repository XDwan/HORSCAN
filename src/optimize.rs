use crate::io::{read_blast_file, save_blast_file, HorScanError, MonAlignmentRow};
#[derive(Clone, Debug)]
struct TypeCount {
    pub mth: i32,
    pub ins: i32,
    pub del: i32,
    pub mis: i32,
}

fn find_useful_target_match(
    blast: &Vec<MonAlignmentRow>,
    current_index: usize,
) -> Result<Vec<usize>, HorScanError> {
    // 初始化一个空的Vec来存储有用的ID
    let mut useful_id = vec![];
    for i in current_index..blast.len() {
        // 确保前方没有占用的mon
        if blast[i].source_mon != "-" {
            break;
        }

        if blast[i].target_mon != "-" && blast[current_index].source_mon == blast[i].target_mon {
            useful_id.push(i);
            break;
        }
    }

    for i in (0..current_index).rev() {
        // 确保后方没有占用的mon
        if blast[i].source_mon != "-" {
            break;
        }
        if blast[i].target_mon != "-" && blast[current_index].source_mon == blast[i].target_mon {
            useful_id.push(i);
            break;
        }
    }
    Ok(useful_id)
}

fn find_useful_source_match(
    blast: &Vec<MonAlignmentRow>,
    current_index: usize,
) -> Result<Vec<usize>, HorScanError> {
    // 初始化一个空的Vec来存储有用的ID
    let mut useful_id = vec![];
    for i in current_index..blast.len() {
        // 确保前方没有占用的mon
        if blast[i].target_mon != "-" {
            break;
        }

        if blast[i].source_mon != "-" && blast[current_index].target_mon == blast[i].source_mon {
            useful_id.push(i);
            break;
        }
    }

    for i in (0..current_index).rev() {
        // 确保后方没有占用的mon
        if blast[i].target_mon != "-" {
            break;
        }
        if blast[i].source_mon != "-" && blast[current_index].target_mon == blast[i].source_mon {
            useful_id.push(i);
            break;
        }
    }
    Ok(useful_id)
}

fn swap_blast_line_source(
    blast: &mut Vec<MonAlignmentRow>,
    i: usize,
    j: usize,
) -> Result<(), HorScanError> {
    // 新建temp row 保存blast[j]的值
    let temp_row = MonAlignmentRow {
        source: blast[j].source.clone(),
        source_start: blast[j].source_start,
        source_end: blast[j].source_end,
        source_mon: blast[j].source_mon.clone(),
        target_start: blast[j].target_start,
        target_end: blast[j].target_end,
        target: blast[j].target.clone(),
        target_mon: blast[j].target_mon.clone(),
        align_type: blast[j].align_type.clone(),
    };
    // 将blast[j]的source的值 替换为 blast[i]的source的值，由于 mon相同 type 替换为MTH
    blast[j].source = blast[i].source.clone();
    blast[j].source_start = blast[i].source_start;
    blast[j].source_end = blast[i].source_end;
    blast[j].source_mon = blast[i].source_mon.clone();
    blast[j].align_type = "MTH".to_string();
    // 将blast[i]的source的值 替换为 temp_row的source的值
    blast[i].source = temp_row.source;
    blast[i].source_start = temp_row.source_start;
    blast[i].source_end = temp_row.source_end;
    blast[i].source_mon = temp_row.source_mon;
    // 如果原来i行的target的mon还存在，则进行标记为DEL， 否则为两个 - - 标记为UNX
    if blast[i].target_mon != "-" {
        blast[i].align_type = "DEL".to_string();
    } else {
        blast[i].align_type = "UNX".to_string();
    }
    Ok(())
}

fn swap_blast_line_target(
    blast: &mut Vec<MonAlignmentRow>,
    i: usize,
    j: usize,
) -> Result<(), HorScanError> {
    // 新建temp row 保存blast[j]的值
    let temp_row = MonAlignmentRow {
        source: blast[j].source.clone(),
        source_start: blast[j].source_start,
        source_end: blast[j].source_end,
        source_mon: blast[j].source_mon.clone(),
        target_start: blast[j].target_start,
        target_end: blast[j].target_end,
        target: blast[j].target.clone(),
        target_mon: blast[j].target_mon.clone(),
        align_type: blast[j].align_type.clone(),
    };
    // 将blast[j]的source的值 替换为 blast[i]的source的值，由于 mon相同 type 替换为MTH
    blast[j].target = blast[i].target.clone();
    blast[j].target_start = blast[i].target_start;
    blast[j].target_end = blast[i].target_end;
    blast[j].target_mon = blast[i].target_mon.clone();
    blast[j].align_type = "MTH".to_string(); // 替换为MTH类型

    // 将blast[i]的source的值 替换为 temp_row的source的值
    blast[i].target = temp_row.target;
    blast[i].target_start = temp_row.target_start;
    blast[i].target_end = temp_row.target_end;
    blast[i].target_mon = temp_row.target_mon.clone();
    // 如果原来i行的target的mon还存在，则进行标记为DEL， 否则为两个 - - 标记为UNX
    if blast[i].target_mon != "-" {
        blast[i].align_type = "INS".to_string();
    } else {
        blast[i].align_type = "UNX".to_string();
    }
    Ok(())
}

// optimize blast file to improve global score
pub fn optimize_blast_file(optimize_path: String) -> Result<(), HorScanError> {
    // 读取blast文件
    let mut blast = read_blast_file(optimize_path.clone())?;

    // 第一遍变量 统计 align_type 的类别数量
    let mut type_count = TypeCount {
        mth: 0,
        ins: 0,
        del: 0,
        mis: 0,
    }; // 初始化结构体
    for row in &blast {
        match row.align_type.as_str() {
            "MTH" => type_count.mth += 1,
            "INS" => type_count.ins += 1,
            "DEL" => type_count.del += 1,
            "MIS" => type_count.mis += 1,
            _ => {} // 忽略未知类型
        }
    }
    println!("TypeCount: {:?}", type_count);

    for i in 0..blast.len() {
        if blast[i].align_type == "MTH" {
            continue;
        }
        // 每一行source 和 target 都能找到匹配并替换，但是不会影响到序列顺序
        let source_swap = find_useful_target_match(&blast, i)?;
        let target_swap = find_useful_source_match(&blast, i)?;

        if source_swap.len() > 0 {
            let j = source_swap[0];
            let _ = swap_blast_line_source(&mut blast, i, j);
        }
        if target_swap.len() > 0 {
            let j = target_swap[0];
            let _ = swap_blast_line_target(&mut blast, i, j);
        }
    }

    // 去除UNX行
    blast.retain(|line| line.align_type != "UNX");
    // 统计optimize 后的结果
    let mut optimize_result = TypeCount {
        mth: 0,
        ins: 0,
        del: 0,
        mis: 0,
    }; 
    for line in &blast {
        match line.align_type.as_str() {
            "MTH" => optimize_result.mth += 1,
            "INS" => optimize_result.ins += 1,
            "DEL" => optimize_result.del += 1,
            "MIS" => optimize_result.mis += 1,
            _ => {}
        }
    }
    println!("optimize : {:?}", optimize_result); // 打印optimize 后的结果

    // 保存修改好的blast文件
    save_blast_file(&blast, optimize_path.clone())?;
    Ok(())
}
