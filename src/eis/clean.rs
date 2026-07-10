use crate::eis::{ImagSignPolicy, ReadOptions, ReadReport, apply_imag_sign, read_spectrum};
use anyhow::{Context, Result, bail};
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};

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

#[derive(Debug, Clone, Serialize)]
pub struct CleanReport {
    pub input_path: PathBuf,
    pub input: ReadReport,
    pub imag_sign: ImagSignPolicy,
    pub drop_positive_imag: bool,
    pub output_point_count: usize,
    pub output_files: Vec<PathBuf>,
}

pub fn clean_file(input: &Path, options: &CleanOptions) -> Result<CleanReport> {
    let stem = input
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("eis");
    let default_root = input
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("result");
    let output_dir = options
        .out_root
        .as_deref()
        .unwrap_or(&default_root)
        .join(stem);
    clean_file_to(input, options, &output_dir)
}

pub fn clean_file_to(
    input: &Path,
    options: &CleanOptions,
    output_dir: &Path,
) -> Result<CleanReport> {
    let mut outcome = read_spectrum(
        input,
        &ReadOptions {
            lenient: options.lenient,
        },
    )?;
    apply_imag_sign(&mut outcome.spectrum, options.imag_sign);
    if options.drop_positive_imag {
        outcome
            .spectrum
            .points
            .retain(|point| point.z_imag_ohm <= 0.0);
    }
    if outcome.spectrum.points.is_empty() {
        bail!("cleaning removed every EIS point from {}", input.display());
    }

    fs::create_dir_all(output_dir)
        .with_context(|| format!("failed to create {}", output_dir.display()))?;
    let csv_path = output_dir.join("cleaned.csv");
    let tsv_path = output_dir.join("cleaned.z60");
    write_points(&csv_path, b',', &outcome.spectrum.points)?;
    write_points(&tsv_path, b'\t', &outcome.spectrum.points)?;
    let report_path = output_dir.join("input_report.json");
    let report = CleanReport {
        input_path: input.to_path_buf(),
        input: outcome.report,
        imag_sign: options.imag_sign,
        drop_positive_imag: options.drop_positive_imag,
        output_point_count: outcome.spectrum.points.len(),
        output_files: vec![csv_path, tsv_path, report_path.clone()],
    };
    fs::write(&report_path, serde_json::to_string_pretty(&report)?)?;
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
