use clap::Parser;

/// Program that processes command line arguments
#[derive(Parser)]
#[command(name = "HORSCAN")]
#[command(version = "1.0")]
#[command(about = "Command line tool for processing genome data", long_about = None)]
pub struct Args {
    /// Input file
    #[arg(short, long,  help = "source bed path")]
    pub source: String,
    #[arg(short, long,  help = "target bed path")]
    pub target: String,
    /// Output file
    #[arg(short, long,  help = "output alignment file prefix")]
    pub output: String,

    #[arg(short, long, default_value = "None", help = "temp file prefix")]
    pub prefix: String,
    /// HORSCAN params
    #[arg(short, long, value_parser, help = "run mode", num_args = 1..)]
    pub mode: Vec<i32>,
}

pub fn parse_args() -> Args {
    Args::parse()
}
