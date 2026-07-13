use anyhow::Result;
use clap::Parser;
use electrochem_tools::batch::{BatchStatus, default_jobs};
use electrochem_tools::eis::{CleanBatchOptions, CleanOptions, ImagSignPolicy, clean_files};
use std::path::PathBuf;

#[derive(Parser, Debug)]
struct Args {
    #[arg(short = 'i', long = "input", required = true, num_args = 1..)]
    input_files: Vec<PathBuf>,
    #[arg(long)]
    lenient: bool,
    #[arg(long, value_enum, default_value_t = ImagSignPolicy::Preserve)]
    imag_sign: ImagSignPolicy,
    #[arg(long)]
    keep_positive_imag: bool,
    #[arg(long)]
    jobs: Option<usize>,
    #[arg(long)]
    fail_fast: bool,
    #[arg(long)]
    overwrite: bool,
    #[arg(
        long,
        help = "Write flat stem-prefixed outputs here instead of beside each input"
    )]
    out_root: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    eprintln!("migration notice: use `eiscli clean -i ...`; clean_eis is a compatibility wrapper");
    let report = clean_files(
        &args.input_files,
        &CleanOptions {
            lenient: args.lenient,
            imag_sign: args.imag_sign,
            drop_positive_imag: !args.keep_positive_imag,
            out_root: args.out_root,
        },
        &CleanBatchOptions {
            jobs: args
                .jobs
                .unwrap_or_else(|| default_jobs(args.input_files.len())),
            fail_fast: args.fail_fast,
            overwrite: args.overwrite,
        },
    )?;
    for item in &report.items {
        match item.status {
            BatchStatus::Success => println!(
                "ok      {} -> {}",
                item.input_path.display(),
                item.output_dir.display()
            ),
            BatchStatus::NotProcessed => println!(
                "not processed {} -> {}",
                item.input_path.display(),
                item.output_dir.display()
            ),
            BatchStatus::Failed => eprintln!(
                "failed  {}: {}",
                item.input_path.display(),
                item.error.as_deref().unwrap_or("unknown error")
            ),
            BatchStatus::Resumed => unreachable!("clean batches cannot resume"),
        }
    }
    println!(
        "batch: {} succeeded, {} failed, {} not processed",
        report.succeeded, report.failed, report.not_processed
    );
    if report.failed > 0 {
        anyhow::bail!("batch completed with {} failed input(s)", report.failed);
    }
    Ok(())
}
