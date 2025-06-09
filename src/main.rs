use io::HorScanError;

use crate::io::Mode;
use crate::horscan::horscan_main; // 引入 horscan_main 函数

mod io;
mod args;
mod horscan;
mod optimize;


fn main() -> Result<(), HorScanError> {
    let params = args::parse_args();
    let mode;
    if params.mode.len() < 3 {
        println!("run with default params");
        mode = Mode {
            // match_method: 0,
            match_score: 4,
            mismatch_score: -5,
            gap_score: -2,
            // jump_max: 10,
        };
        println!("{:?}", mode);
    } else {
        mode = Mode {
            // match_method: params.mode[0],
            match_score: params.mode[0],
            mismatch_score: params.mode[1] * -1,
            gap_score: params.mode[2] * -1,
            // jump_max: params.mode[4] as usize,
        };
        println!("{:?}", mode);
    }
    horscan_main(params.source, params.target, params.output, mode)?;
    Ok(())
}