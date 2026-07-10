use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

#[derive(Debug, Clone)]
pub struct BatchOptions {
    pub jobs: usize,
    pub fail_fast: bool,
    pub overwrite: bool,
    pub resume: bool,
    pub out_root: PathBuf,
    pub output_suffix: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BatchStatus {
    Success,
    Failed,
    Resumed,
    NotProcessed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchItem {
    pub input_index: usize,
    pub input_path: PathBuf,
    pub output_dir: PathBuf,
    pub status: BatchStatus,
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchReport {
    pub items: Vec<BatchItem>,
    pub succeeded: usize,
    pub failed: usize,
    pub skipped: usize,
    pub resumed: usize,
    pub not_processed: usize,
}

pub fn default_jobs(input_count: usize) -> usize {
    std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
        .min(input_count.max(1))
}

pub fn run_batch<F>(inputs: &[PathBuf], options: &BatchOptions, processor: F) -> Result<BatchReport>
where
    F: Fn(&Path, &Path) -> Result<()> + Sync,
{
    run_batch_with_resume(
        inputs,
        options,
        |_, _| bail!("resume validation is not configured for this batch"),
        processor,
    )
}

pub fn output_dir_for(input: &Path, options: &BatchOptions) -> Result<PathBuf> {
    let stem = input
        .file_stem()
        .and_then(|value| value.to_str())
        .ok_or_else(|| {
            anyhow::anyhow!(
                "input path has no valid UTF-8 file stem: {}",
                input.display()
            )
        })?;

    let logical_stem = stem
        .strip_suffix("_cleaned")
        .or_else(|| stem.strip_suffix("-cleaned"))
        .unwrap_or(stem);

    Ok(options
        .out_root
        .join(format!("{logical_stem}_{}", options.output_suffix)))
}

pub fn run_batch_with_resume<R, F>(
    inputs: &[PathBuf],
    options: &BatchOptions,
    resume_validator: R,
    processor: F,
) -> Result<BatchReport>
where
    R: Fn(&Path, &Path) -> Result<bool> + Sync,
    F: Fn(&Path, &Path) -> Result<()> + Sync,
{
    if inputs.is_empty() {
        bail!("batch requires at least one input file");
    }
    if options.jobs == 0 {
        bail!("--jobs must be at least 1");
    }
    if options.resume && options.overwrite {
        bail!("--resume and --overwrite cannot be used together");
    }
    fs::create_dir_all(&options.out_root)?;
    let next = AtomicUsize::new(0);
    let stop = AtomicBool::new(false);
    let results: Mutex<Vec<Option<BatchItem>>> = Mutex::new(vec![None; inputs.len()]);
    let jobs = options.jobs.min(inputs.len());

    let output_dirs = inputs
        .iter()
        .map(|input| output_dir_for(input, options))
        .collect::<Result<Vec<_>>>()?;

    std::thread::scope(|scope| {
        for _ in 0..jobs {
            scope.spawn(|| {
                loop {
                    if options.fail_fast && stop.load(Ordering::Acquire) {
                        break;
                    }
                    let index = next.fetch_add(1, Ordering::AcqRel);
                    if index >= inputs.len() {
                        break;
                    }
                    let input = &inputs[index];
                    let output_dir = output_dirs[index].clone();
                    let result = prepare_output(input, &output_dir, options, &resume_validator)
                        .and_then(|should_run| {
                            should_run
                                .then(|| processor(input, &output_dir))
                                .transpose()
                        });
                    let item = match result {
                        Ok(Some(())) => BatchItem {
                            input_index: index,
                            input_path: input.clone(),
                            output_dir,
                            status: BatchStatus::Success,
                            error: None,
                        },
                        Ok(None) => BatchItem {
                            input_index: index,
                            input_path: input.clone(),
                            output_dir,
                            status: BatchStatus::Resumed,
                            error: None,
                        },
                        Err(error) => {
                            stop.store(true, Ordering::Release);
                            BatchItem {
                                input_index: index,
                                input_path: input.clone(),
                                output_dir,
                                status: BatchStatus::Failed,
                                error: Some(format!("{error:#}")),
                            }
                        }
                    };
                    results.lock().expect("batch result mutex poisoned")[index] = Some(item);
                }
            });
        }
    });

    let mut locked = results.into_inner().expect("batch result mutex poisoned");
    let items = inputs
        .iter()
        .enumerate()
        .map(|(index, input)| {
            locked[index].take().unwrap_or_else(|| BatchItem {
                input_index: index,
                input_path: input.clone(),
                output_dir: options.out_root.join(output_dirs[index].file_name().unwrap()),
                status: BatchStatus::NotProcessed,
                error: Some("not processed because --fail-fast stopped the batch".to_string()),
            })
        })
        .collect::<Vec<_>>();
    let report = BatchReport {
        succeeded: items
            .iter()
            .filter(|item| item.status == BatchStatus::Success)
            .count(),
        failed: items
            .iter()
            .filter(|item| item.status == BatchStatus::Failed)
            .count(),
        skipped: items
            .iter()
            .filter(|item| {
                matches!(
                    item.status,
                    BatchStatus::Resumed | BatchStatus::NotProcessed
                )
            })
            .count(),
        resumed: items
            .iter()
            .filter(|item| item.status == BatchStatus::Resumed)
            .count(),
        not_processed: items
            .iter()
            .filter(|item| item.status == BatchStatus::NotProcessed)
            .count(),
        items,
    };
    write_batch_summary(&options.out_root.join("batch_summary.csv"), &report)?;
    Ok(report)
}

fn prepare_output<R>(
    input: &Path,
    output_dir: &Path,
    options: &BatchOptions,
    resume_validator: &R,
) -> Result<bool>
where
    R: Fn(&Path, &Path) -> Result<bool>,
{
    if output_dir.exists() {
        if options.resume {
            return resume_validator(input, output_dir);
        }
        if !options.overwrite {
            bail!(
                "output directory {} already exists; use --overwrite or choose another --out-root",
                output_dir.display()
            );
        }
        fs::remove_dir_all(output_dir)?;
    }
    fs::create_dir_all(output_dir)?;
    Ok(true)
}

fn write_batch_summary(path: &Path, report: &BatchReport) -> Result<()> {
    let mut writer = csv::Writer::from_path(path)?;
    for item in &report.items {
        writer.serialize(item)?;
    }
    writer.flush()?;
    Ok(())
}
