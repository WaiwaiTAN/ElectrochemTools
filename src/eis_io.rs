use crate::types::EisData;
use anyhow::{Context, Result, bail};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Copy)]
struct ColumnMap {
    frequency: usize,
    z_real: usize,
    z_imag: usize,
}

pub fn read_eis(path: &Path) -> Result<EisData> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read EIS file {}", path.display()))?;
    let mut rows = Vec::new();

    let parse_plan = detect_parse_plan(&content)
        .with_context(|| format!("could not detect EIS columns in {}", path.display()))?;

    for (line_idx, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || parse_plan.header_line.is_some_and(|idx| line_idx <= idx) {
            continue;
        }

        let fields = split_line(trimmed, parse_plan.delimiter);
        let Some(row) = parse_row(&fields, parse_plan.column_map) else {
            if parse_plan.has_header {
                continue;
            }
            bail!("failed to parse numeric EIS row {}: {}", line_idx + 1, line);
        };
        if row.0.is_finite() && row.0 > 0.0 && row.1.is_finite() && row.2.is_finite() {
            rows.push(row);
        }
    }

    if rows.is_empty() {
        bail!("no valid EIS rows found in {}", path.display());
    }

    rows.sort_by(|a, b| b.0.total_cmp(&a.0));
    EisData::new(
        rows.iter().map(|row| row.0).collect(),
        rows.iter().map(|row| row.1).collect(),
        rows.iter().map(|row| row.2).collect(),
    )
}

#[derive(Debug, Clone, Copy)]
struct ParsePlan {
    delimiter: char,
    column_map: ColumnMap,
    header_line: Option<usize>,
    has_header: bool,
}

fn detect_parse_plan(content: &str) -> Result<ParsePlan> {
    for (line_idx, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let delimiter = detect_delimiter(trimmed);
        let fields = split_line(trimmed, delimiter);
        if fields_are_numeric(&fields) {
            continue;
        }
        if let Ok(column_map) = detect_columns(&fields) {
            return Ok(ParsePlan {
                delimiter,
                column_map,
                header_line: Some(line_idx),
                has_header: true,
            });
        }
    }

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let delimiter = detect_delimiter(trimmed);
        let fields = split_line(trimmed, delimiter);
        if fields_are_numeric(&fields) {
            return Ok(ParsePlan {
                delimiter,
                column_map: detect_numeric_columns(&fields)?,
                header_line: None,
                has_header: false,
            });
        }
    }

    bail!("no EIS header or numeric data row found")
}

pub fn read_eis_with_cleaning(path: &Path, drop_positive_imag: bool) -> Result<EisData> {
    let mut data = read_eis(path)?;
    if drop_positive_imag {
        let rows: Vec<(f64, f64, f64)> = data
            .frequency_hz
            .iter()
            .copied()
            .zip(data.z_real.iter().copied())
            .zip(data.z_imag.iter().copied())
            .map(|((f, re), im)| (f, re, im))
            .filter(|(_, _, im)| *im <= 0.0)
            .collect();
        if rows.is_empty() {
            bail!("all EIS rows were removed by positive-imaginary filtering");
        }
        data = EisData::new(
            rows.iter().map(|row| row.0).collect(),
            rows.iter().map(|row| row.1).collect(),
            rows.iter().map(|row| row.2).collect(),
        )?;
    } else {
        data.sort_by_frequency_desc();
    }
    Ok(data)
}

pub fn write_impedance_csv(
    path: &Path,
    data: &EisData,
    z_fit: &[num_complex::Complex<f64>],
) -> Result<()> {
    if data.len() != z_fit.len() {
        bail!("cannot write impedance CSV: fit length does not match data length");
    }
    let mut writer = csv::Writer::from_path(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    writer.write_record([
        "frequency",
        "Z_real_exp",
        "Z_imag_exp",
        "Z_real_fit",
        "Z_imag_fit",
        "residual_real",
        "residual_imag",
    ])?;
    for (((&freq, &re), &im), fit) in data
        .frequency_hz
        .iter()
        .zip(&data.z_real)
        .zip(&data.z_imag)
        .zip(z_fit)
    {
        writer.serialize((freq, re, im, fit.re, fit.im, re - fit.re, im - fit.im))?;
    }
    writer.flush()?;
    Ok(())
}

fn detect_delimiter(line: &str) -> char {
    let comma = line.matches(',').count();
    let tab = line.matches('\t').count();
    let semicolon = line.matches(';').count();
    if tab >= comma && tab >= semicolon && tab > 0 {
        '\t'
    } else if semicolon > comma {
        ';'
    } else {
        ','
    }
}

fn split_line(line: &str, delimiter: char) -> Vec<String> {
    line.split(delimiter)
        .map(|field| field.trim().trim_matches('"').to_string())
        .collect()
}

fn fields_are_numeric(fields: &[String]) -> bool {
    fields.iter().any(|field| field.parse::<f64>().is_ok())
        && fields
            .iter()
            .filter(|field| !field.trim().is_empty())
            .all(|field| field.parse::<f64>().is_ok())
}

fn detect_columns(headers: &[String]) -> Result<ColumnMap> {
    let normalized: Vec<String> = headers
        .iter()
        .map(|header| normalize_header(header))
        .collect();
    let frequency = find_first(&normalized, &["freq", "frequency", "frequencyhz", "freqhz"])
        .or_else(|| normalized.iter().position(|h| h.contains("freq")))
        .ok_or_else(|| anyhow::anyhow!("missing frequency column"))?;
    let z_real = find_first(&normalized, &["zreal", "zre", "rez", "zr", "zprime", "z"])
        .or_else(|| headers.iter().position(|h| is_real_impedance_header(h)))
        .or_else(|| normalized.iter().position(|h| h.contains("real")))
        .ok_or_else(|| anyhow::anyhow!("missing real impedance column"))?;
    let z_imag = find_first(
        &normalized,
        &["zimag", "zim", "imz", "zi", "zdoubleprime", "izr"],
    )
    .or_else(|| headers.iter().position(|h| is_imag_impedance_header(h)))
    .or_else(|| normalized.iter().position(|h| h.contains("imag")))
    .ok_or_else(|| anyhow::anyhow!("missing imaginary impedance column"))?;
    Ok(ColumnMap {
        frequency,
        z_real,
        z_imag,
    })
}

fn is_real_impedance_header(header: &str) -> bool {
    let value = header.trim().to_ascii_lowercase();
    value == "z'" || value.starts_with("z'(") || value == "zr" || value == "z_real"
}

fn is_imag_impedance_header(header: &str) -> bool {
    let value = header.trim().to_ascii_lowercase();
    value == "z''" || value.starts_with("z''(") || value == "izr" || value == "z_imag"
}

fn detect_numeric_columns(fields: &[String]) -> Result<ColumnMap> {
    let numeric: Vec<usize> = fields
        .iter()
        .enumerate()
        .filter_map(|(idx, field)| field.parse::<f64>().ok().map(|_| idx))
        .collect();
    if numeric.len() < 3 {
        bail!("no-header EIS data must contain at least three numeric columns");
    }
    Ok(ColumnMap {
        frequency: numeric[0],
        z_real: numeric[1],
        z_imag: numeric[2],
    })
}

fn parse_row(fields: &[String], map: ColumnMap) -> Option<(f64, f64, f64)> {
    Some((
        fields.get(map.frequency)?.parse().ok()?,
        fields.get(map.z_real)?.parse().ok()?,
        fields.get(map.z_imag)?.parse().ok()?,
    ))
}

fn normalize_header(header: &str) -> String {
    header
        .to_ascii_lowercase()
        .chars()
        .filter(|ch| ch.is_ascii_alphanumeric())
        .collect()
}

fn find_first(headers: &[String], candidates: &[&str]) -> Option<usize> {
    candidates
        .iter()
        .find_map(|candidate| headers.iter().position(|header| header == candidate))
}
