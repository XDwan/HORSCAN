mod args;
mod errors;
mod types;
mod io;
mod score;
mod align;
mod horscan;

use anyhow::Result;

fn main() -> Result<()> {
    // init logger: RUST_LOG=info ./horscanv ...
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cfg = args::Args::parse();
    horscan::horscan_main(cfg)
}
