use anyhow::Result;
use clap::Parser;
use electrochem_tools::eis::{CleanOptions, clean_file};
use std::path::PathBuf;

#[derive(Parser, Debug)]
struct Args {
    #[arg(short = 'i', long = "input", required = true, num_args = 1..)]
    input_files: Vec<PathBuf>,
    #[arg(long)]
    lenient: bool,
    #[arg(long)]
    out_root: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    eprintln!("migration notice: use `eiscli clean -i ...`; clean_eis is a compatibility wrapper");
    for input in &args.input_files {
        let report = clean_file(
            input,
            &CleanOptions {
                lenient: args.lenient,
                out_root: args.out_root.clone(),
                ..CleanOptions::default()
            },
        )?;
        println!(
            "cleaned {}: read {}, skipped {}, wrote {} points",
            input.display(),
            report.valid_row_count,
            report.skipped_row_count,
            report.output_point_count
        );
    }
    Ok(())
}
