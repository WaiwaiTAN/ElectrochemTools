use crate::batch::{BatchItem, BatchReport, BatchStatus};
use crate::eis::{EisFormat, ImagSignPolicy, ReadOptions, apply_imag_sign, read_spectrum};
use anyhow::{Context, Result, bail};
use serde::Serialize;
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

#[derive(Debug, Clone)]
pub struct CleanOptions {
    pub lenient: bool,
    pub imag_sign: ImagSignPolicy,
    pub drop_positive_imag: bool,
    pub out_root: Option<PathBuf>,
}

impl Default for CleanOptions {
    fn default() -> Self {
        Self {
            lenient: false,
            imag_sign: ImagSignPolicy::Preserve,
            drop_positive_imag: true,
            out_root: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CleanBatchOptions {
    pub jobs: usize,
    pub fail_fast: bool,
    pub overwrite: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct CleanReport {
    pub input_path: PathBuf,
    pub detected_input_format: EisFormat,
    pub mode: String,
    pub imag_sign: ImagSignPolicy,
    pub drop_positive_imag: bool,
    pub original_row_count: usize,
    pub valid_row_count: usize,
    pub skipped_row_count: usize,
    pub skipped_by_reason: std::collections::BTreeMap<String, usize>,
    pub removed_positive_imag_count: usize,
    pub output_point_count: usize,
    pub output_files: Vec<PathBuf>,
}

pub fn clean_file(input: &Path, options: &CleanOptions) -> Result<CleanReport> {
    let output_files = clean_output_files(input, options);
    clean_file_to_paths(
        input,
        options,
        &output_files[0],
        &output_files[1],
        &output_files[2],
    )
}

pub fn clean_files(
    inputs: &[PathBuf],
    options: &CleanOptions,
    batch: &CleanBatchOptions,
) -> Result<BatchReport> {
    if inputs.is_empty() {
        bail!("cleaning requires at least one input file");
    }
    if batch.jobs == 0 {
        bail!("--jobs must be at least 1");
    }

    let output_files = inputs
        .iter()
        .map(|input| clean_output_files(input, options))
        .collect::<Vec<_>>();
    let mut unique_outputs = HashSet::new();
    for (input, files) in inputs.iter().zip(&output_files) {
        for path in files {
            if !unique_outputs.insert(path.clone()) {
                bail!(
                    "clean inputs map to the same output path {}; use distinct input stems or output locations (conflict at {})",
                    path.display(),
                    input.display()
                );
            }
        }
    }

    let next = AtomicUsize::new(0);
    let stop = AtomicBool::new(false);
    let results: Mutex<Vec<Option<BatchItem>>> = Mutex::new(vec![None; inputs.len()]);
    let jobs = batch.jobs.min(inputs.len());
    std::thread::scope(|scope| {
        for _ in 0..jobs {
            scope.spawn(|| {
                loop {
                    if batch.fail_fast && stop.load(Ordering::Acquire) {
                        break;
                    }
                    let index = next.fetch_add(1, Ordering::AcqRel);
                    if index >= inputs.len() {
                        break;
                    }
                    let input = &inputs[index];
                    let files = &output_files[index];
                    let output_dir = files[0]
                        .parent()
                        .unwrap_or_else(|| Path::new("."))
                        .to_path_buf();
                    let result = if !batch.overwrite && files.iter().any(|path| path.exists()) {
                        Err(anyhow::anyhow!(
                            "clean output for {} already exists; use --overwrite",
                            input.display()
                        ))
                    } else {
                        clean_file(input, options).map(|_| ())
                    };
                    let item = match result {
                        Ok(()) => BatchItem {
                            input_index: index,
                            input_path: input.clone(),
                            output_dir,
                            status: BatchStatus::Success,
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
                    results.lock().expect("clean result mutex poisoned")[index] = Some(item);
                }
            });
        }
    });

    let mut locked = results.into_inner().expect("clean result mutex poisoned");
    let items = inputs
        .iter()
        .enumerate()
        .map(|(index, input)| {
            locked[index].take().unwrap_or_else(|| BatchItem {
                input_index: index,
                input_path: input.clone(),
                output_dir: output_files[index][0]
                    .parent()
                    .unwrap_or_else(|| Path::new("."))
                    .to_path_buf(),
                status: BatchStatus::NotProcessed,
                error: Some("not processed because --fail-fast stopped the batch".to_string()),
            })
        })
        .collect::<Vec<_>>();
    Ok(BatchReport {
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
            .filter(|item| item.status == BatchStatus::NotProcessed)
            .count(),
        resumed: 0,
        not_processed: items
            .iter()
            .filter(|item| item.status == BatchStatus::NotProcessed)
            .count(),
        items,
    })
}

fn clean_output_files(input: &Path, options: &CleanOptions) -> [PathBuf; 3] {
    let stem = input
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("eis");
    let default_root = input
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let output_dir = options.out_root.as_deref().unwrap_or(default_root);
    [
        output_dir.join(format!("{stem}_cleaned.csv")),
        output_dir.join(format!("{stem}_cleaned.z60")),
        output_dir.join(format!("{stem}_clean_state.json")),
    ]
}

pub fn clean_file_to(
    input: &Path,
    options: &CleanOptions,
    output_dir: &Path,
) -> Result<CleanReport> {
    clean_file_to_paths(
        input,
        options,
        &output_dir.join("cleaned.csv"),
        &output_dir.join("cleaned.z60"),
        &output_dir.join("input_report.json"),
    )
}

fn clean_file_to_paths(
    input: &Path,
    options: &CleanOptions,
    csv_path: &Path,
    tsv_path: &Path,
    report_path: &Path,
) -> Result<CleanReport> {
    let mut outcome = read_spectrum(
        input,
        &ReadOptions {
            lenient: options.lenient,
        },
    )?;
    let detected_input_format = outcome.spectrum.metadata.source_format;
    let valid_row_count = outcome.spectrum.points.len();
    apply_imag_sign(&mut outcome.spectrum, options.imag_sign);
    let before_positive_filter = outcome.spectrum.points.len();
    if options.drop_positive_imag {
        outcome
            .spectrum
            .points
            .retain(|point| point.z_imag_ohm <= 0.0);
    }
    if outcome.spectrum.points.is_empty() {
        bail!("cleaning removed every EIS point from {}", input.display());
    }

    let output_dir = csv_path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_dir)
        .with_context(|| format!("failed to create {}", output_dir.display()))?;
    write_points(csv_path, b',', &outcome.spectrum.points)?;
    write_points(tsv_path, b'\t', &outcome.spectrum.points)?;
    let report = CleanReport {
        input_path: input.to_path_buf(),
        detected_input_format,
        mode: if options.lenient { "lenient" } else { "strict" }.to_string(),
        imag_sign: options.imag_sign,
        drop_positive_imag: options.drop_positive_imag,
        original_row_count: outcome.report.total_rows,
        valid_row_count,
        skipped_row_count: outcome.report.rows_skipped,
        skipped_by_reason: outcome.report.skipped_by_reason,
        removed_positive_imag_count: before_positive_filter - outcome.spectrum.points.len(),
        output_point_count: outcome.spectrum.points.len(),
        output_files: vec![
            csv_path.to_path_buf(),
            tsv_path.to_path_buf(),
            report_path.to_path_buf(),
        ],
    };
    fs::write(report_path, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
}

fn write_points(path: &Path, delimiter: u8, points: &[crate::eis::EisPoint]) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_path(path)?;
    writer.write_record(["frequency_hz", "z_real_ohm", "z_imag_ohm"])?;
    for point in points {
        writer.serialize((point.frequency_hz, point.z_real_ohm, point.z_imag_ohm))?;
    }
    writer.flush()?;
    Ok(())
}
