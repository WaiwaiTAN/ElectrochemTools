use crate::eis::data::{EisMetadata, EisPoint, EisSpectrum};
use crate::eis::format::EisFormat;
use crate::eis::validation::{detect_imaginary_convention, validate_point};
use anyhow::{Context, Result, anyhow, bail};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Copy, Default)]
pub struct ReadOptions {
    pub lenient: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReadReport {
    pub total_rows: usize,
    pub rows_read: usize,
    pub rows_skipped: usize,
    pub skipped_by_reason: BTreeMap<String, usize>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct ReadOutcome {
    pub spectrum: EisSpectrum,
    pub report: ReadReport,
}

#[derive(Debug, Clone, Copy)]
struct ColumnMap {
    frequency: usize,
    z_real: usize,
    z_imag: usize,
}

#[derive(Debug, Clone, Copy)]
struct ParsePlan {
    delimiter: char,
    columns: ColumnMap,
    first_data_line: usize,
    has_header: bool,
}

pub fn read_spectrum(path: &Path, options: &ReadOptions) -> Result<ReadOutcome> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read EIS file {}", path.display()))?;
    let plan = detect_parse_plan(&content)
        .with_context(|| format!("could not detect EIS columns in {}", path.display()))?;
    let mut points = Vec::new();
    let mut skipped = BTreeMap::new();
    let mut total_rows = 0;
    let mut seen_frequency = HashSet::new();

    for (line_index, line) in content.lines().enumerate().skip(plan.first_data_line) {
        if line.trim().is_empty() || line.trim().eq_ignore_ascii_case("End Comments") {
            continue;
        }
        total_rows += 1;
        let fields = split_line(line.trim(), plan.delimiter);
        let point = match parse_point(&fields, plan.columns) {
            Ok(point) => point,
            Err(error) => {
                handle_bad_row(options, &mut skipped, "parse-error", line_index, error)?;
                continue;
            }
        };
        if let Err(error) = validate_point(&point) {
            let reason = if !point.frequency_hz.is_finite()
                || !point.z_real_ohm.is_finite()
                || !point.z_imag_ohm.is_finite()
            {
                "non-finite-value"
            } else {
                "non-positive-frequency"
            };
            handle_bad_row(options, &mut skipped, reason, line_index, error)?;
            continue;
        }
        if !seen_frequency.insert(point.frequency_hz.to_bits()) {
            handle_bad_row(
                options,
                &mut skipped,
                "duplicate-frequency",
                line_index,
                anyhow!("duplicate frequency {} Hz", point.frequency_hz),
            )?;
            continue;
        }
        points.push(point);
    }
    if points.is_empty() {
        bail!("no valid EIS rows found in {}", path.display());
    }

    let non_monotonic = points
        .windows(2)
        .any(|pair| pair[0].frequency_hz <= pair[1].frequency_hz);
    let mut warnings = Vec::new();
    if non_monotonic {
        warnings.push("input frequency is not strictly descending; points were sorted".to_string());
    }
    points.sort_by(|left, right| right.frequency_hz.total_cmp(&left.frequency_hz));
    let convention = detect_imaginary_convention(&points);
    if convention == crate::eis::data::ImaginaryConvention::PositiveCapacitive {
        warnings
            .push("most imaginary impedance values are positive; sign was preserved".to_string());
    }
    let rows_skipped = skipped.values().sum();
    let source_format = detect_format(path, plan);
    let metadata = EisMetadata {
        source_path: Some(path.to_path_buf()),
        source_format,
        imaginary_convention: convention,
        original_point_count: total_rows,
        preprocessing: Vec::new(),
        warnings: warnings.clone(),
    };
    Ok(ReadOutcome {
        report: ReadReport {
            total_rows,
            rows_read: points.len(),
            rows_skipped,
            skipped_by_reason: skipped,
            warnings,
        },
        spectrum: EisSpectrum { points, metadata },
    })
}

fn handle_bad_row(
    options: &ReadOptions,
    skipped: &mut BTreeMap<String, usize>,
    reason: &str,
    zero_based_line: usize,
    error: anyhow::Error,
) -> Result<()> {
    if !options.lenient {
        bail!(
            "invalid EIS row {} ({reason}): {error}",
            zero_based_line + 1
        );
    }
    *skipped.entry(reason.to_string()).or_default() += 1;
    Ok(())
}

fn detect_parse_plan(content: &str) -> Result<ParsePlan> {
    for (line_index, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let delimiter = detect_delimiter(trimmed);
        let fields = split_line(trimmed, delimiter);
        if fields.iter().all(|field| field.parse::<f64>().is_ok()) {
            if fields.len() < 3 {
                bail!("headerless EIS rows require at least three columns");
            }
            return Ok(ParsePlan {
                delimiter,
                columns: ColumnMap {
                    frequency: 0,
                    z_real: 1,
                    z_imag: 2,
                },
                first_data_line: line_index,
                has_header: false,
            });
        }
        match detect_columns(&fields) {
            Ok(columns) => {
                return Ok(ParsePlan {
                    delimiter,
                    columns,
                    first_data_line: line_index + 1,
                    has_header: true,
                });
            }
            Err(error) if looks_like_eis_header(&fields) => return Err(error),
            Err(_) => {}
        }
    }
    bail!("no EIS header or numeric data row found")
}

fn looks_like_eis_header(fields: &[String]) -> bool {
    let normalized: Vec<String> = fields.iter().map(|field| normalize_header(field)).collect();
    let has_frequency = normalized.iter().any(|field| field.contains("freq"));
    let has_real = normalized
        .iter()
        .any(|field| field.contains("zreal") || matches!(field.as_str(), "zr" | "zprime"));
    let has_imaginary = normalized
        .iter()
        .any(|field| field.contains("zimag") || matches!(field.as_str(), "izr" | "zdoubleprime"));
    [has_frequency, has_real, has_imaginary]
        .into_iter()
        .filter(|present| *present)
        .count()
        >= 2
}

fn detect_columns(headers: &[String]) -> Result<ColumnMap> {
    let normalized: Vec<String> = headers
        .iter()
        .map(|value| normalize_header(value))
        .collect();
    Ok(ColumnMap {
        frequency: unique_column(
            &normalized,
            |value| {
                matches!(value, "freq" | "frequency" | "frequencyhz" | "freqhz")
                    || value.contains("freq")
            },
            "frequency",
        )?,
        z_real: unique_indexes(
            headers
                .iter()
                .zip(&normalized)
                .enumerate()
                .filter_map(|(index, (raw, value))| {
                    let raw = raw.trim().to_ascii_lowercase();
                    (matches!(
                        value.as_str(),
                        "zreal" | "zre" | "rez" | "zr" | "zprime" | "z"
                    ) || value.contains("real")
                        || raw == "z'"
                        || (raw.starts_with("z'(") && !raw.starts_with("z''")))
                    .then_some(index)
                })
                .collect(),
            "real impedance",
        )?,
        z_imag: unique_indexes(
            headers
                .iter()
                .zip(&normalized)
                .enumerate()
                .filter_map(|(index, (raw, value))| {
                    let raw = raw.trim().to_ascii_lowercase();
                    (matches!(
                        value.as_str(),
                        "zimag" | "zim" | "imz" | "zi" | "zdoubleprime" | "izr"
                    ) || value.contains("imag")
                        || raw == "z''"
                        || raw.starts_with("z''("))
                    .then_some(index)
                })
                .collect(),
            "imaginary impedance",
        )?,
    })
}

fn unique_column(
    headers: &[String],
    predicate: impl Fn(&str) -> bool,
    label: &str,
) -> Result<usize> {
    let matches: Vec<usize> = headers
        .iter()
        .enumerate()
        .filter_map(|(index, header)| predicate(header).then_some(index))
        .collect();
    unique_indexes(matches, label)
}

fn unique_indexes(matches: Vec<usize>, label: &str) -> Result<usize> {
    match matches.as_slice() {
        [index] => Ok(*index),
        [] => bail!("missing {label} column"),
        _ => bail!("ambiguous {label} columns at indexes {matches:?}"),
    }
}

fn parse_point(fields: &[String], map: ColumnMap) -> Result<EisPoint> {
    let parse = |index: usize, label: &str| -> Result<f64> {
        fields
            .get(index)
            .ok_or_else(|| anyhow!("missing {label} field"))?
            .parse::<f64>()
            .with_context(|| format!("cannot parse {label}"))
    };
    Ok(EisPoint {
        frequency_hz: parse(map.frequency, "frequency")?,
        z_real_ohm: parse(map.z_real, "real impedance")?,
        z_imag_ohm: parse(map.z_imag, "imaginary impedance")?,
    })
}

fn detect_format(path: &Path, plan: ParsePlan) -> EisFormat {
    let extension = path
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    match (extension.as_str(), plan.has_header, plan.delimiter) {
        ("z60", _, _) => EisFormat::CorrTestZ60,
        ("txt", true, _) => EisFormat::CorrTestText,
        ("csv", false, ',') if plan.columns.frequency == 0 => EisFormat::LegacyCleanCsv,
        (_, false, _) => EisFormat::HeaderlessThreeColumn,
        (_, _, '\t') => EisFormat::Tsv,
        _ => EisFormat::Csv,
    }
}

fn detect_delimiter(line: &str) -> char {
    let counts = [
        (',', line.matches(',').count()),
        ('\t', line.matches('\t').count()),
        (';', line.matches(';').count()),
    ];
    counts
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .filter(|(_, count)| *count > 0)
        .map(|(delimiter, _)| delimiter)
        .unwrap_or(',')
}

fn split_line(line: &str, delimiter: char) -> Vec<String> {
    line.split(delimiter)
        .map(|field| field.trim().trim_matches('"').to_string())
        .collect()
}

fn normalize_header(header: &str) -> String {
    header
        .to_ascii_lowercase()
        .chars()
        .filter(|character| character.is_ascii_alphanumeric())
        .collect()
}
