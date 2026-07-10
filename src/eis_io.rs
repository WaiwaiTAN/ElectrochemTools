//! Compatibility API for the original vector-based EIS interface.
//!
//! New code should use [`crate::eis`] so metadata and input diagnostics are retained.

use crate::eis::{ReadOptions, read_spectrum};
use crate::types::EisData;
use anyhow::{Context, Result, bail};
use std::path::Path;

pub fn read_eis(path: &Path) -> Result<EisData> {
    let outcome = read_spectrum(path, &ReadOptions::default())?;
    spectrum_to_legacy(outcome.spectrum)
}

pub fn read_eis_with_cleaning(path: &Path, drop_positive_imag: bool) -> Result<EisData> {
    let outcome = read_spectrum(path, &ReadOptions::default())?;
    let mut spectrum = outcome.spectrum;
    if drop_positive_imag {
        spectrum.points.retain(|point| point.z_imag_ohm <= 0.0);
        if spectrum.points.is_empty() {
            bail!("all EIS rows were removed by positive-imaginary filtering");
        }
    }
    spectrum_to_legacy(spectrum)
}

fn spectrum_to_legacy(spectrum: crate::eis::EisSpectrum) -> Result<EisData> {
    EisData::new(
        spectrum
            .points
            .iter()
            .map(|point| point.frequency_hz)
            .collect(),
        spectrum
            .points
            .iter()
            .map(|point| point.z_real_ohm)
            .collect(),
        spectrum
            .points
            .iter()
            .map(|point| point.z_imag_ohm)
            .collect(),
    )
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
