use anyhow::{Result, bail};
use num_complex::Complex;
use serde::Serialize;

#[derive(Debug, Clone)]
pub struct EisData {
    pub frequency_hz: Vec<f64>,
    pub z_real: Vec<f64>,
    pub z_imag: Vec<f64>,
}

impl EisData {
    pub fn new(frequency_hz: Vec<f64>, z_real: Vec<f64>, z_imag: Vec<f64>) -> Result<Self> {
        if frequency_hz.len() != z_real.len() || frequency_hz.len() != z_imag.len() {
            bail!("EIS vectors must have equal length");
        }
        if frequency_hz.is_empty() {
            bail!("EIS data is empty");
        }
        Ok(Self {
            frequency_hz,
            z_real,
            z_imag,
        })
    }

    pub fn len(&self) -> usize {
        self.frequency_hz.len()
    }

    pub fn is_empty(&self) -> bool {
        self.frequency_hz.is_empty()
    }

    pub fn complex_impedance(&self) -> Vec<Complex<f64>> {
        self.z_real
            .iter()
            .zip(&self.z_imag)
            .map(|(&re, &im)| Complex::new(re, im))
            .collect()
    }

    pub fn sort_by_frequency_desc(&mut self) {
        let mut rows: Vec<(f64, f64, f64)> = self
            .frequency_hz
            .iter()
            .copied()
            .zip(self.z_real.iter().copied())
            .zip(self.z_imag.iter().copied())
            .map(|((f, re), im)| (f, re, im))
            .collect();
        rows.sort_by(|a, b| b.0.total_cmp(&a.0));
        self.frequency_hz = rows.iter().map(|row| row.0).collect();
        self.z_real = rows.iter().map(|row| row.1).collect();
        self.z_imag = rows.iter().map(|row| row.2).collect();
    }

    pub fn flip_imag(&mut self) {
        for value in &mut self.z_imag {
            *value = -*value;
        }
    }

    pub fn warn_if_imag_mostly_positive(&self) {
        let positives = self.z_imag.iter().filter(|&&value| value > 0.0).count();
        if positives * 2 > self.z_imag.len() {
            eprintln!(
                "warning: most Z_imag values are positive; sign is preserved. Use --flip-imag to invert explicitly."
            );
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct FitMetrics {
    pub rmse_real: f64,
    pub rmse_imag: f64,
    pub rmse_complex: f64,
    pub rmse_magnitude: f64,
    pub relative_rmse: f64,
    pub relative_rmse_percent: f64,
}

pub fn calculate_fit_metrics(exp: &EisData, z_fit: &[Complex<f64>]) -> Result<FitMetrics> {
    if exp.len() != z_fit.len() {
        bail!("fit vector length does not match EIS data length");
    }
    let n = exp.len() as f64;
    let mut ss_re = 0.0;
    let mut ss_im = 0.0;
    let mut ss_complex = 0.0;
    let mut ss_reference = 0.0;
    for ((&re, &im), fit) in exp.z_real.iter().zip(&exp.z_imag).zip(z_fit) {
        let dr = re - fit.re;
        let di = im - fit.im;
        ss_re += dr * dr;
        ss_im += di * di;
        ss_complex += dr * dr + di * di;
        ss_reference += re * re + im * im;
    }
    let rmse_real = (ss_re / n).sqrt();
    let rmse_imag = (ss_im / n).sqrt();
    let rmse_magnitude = (ss_complex / n).sqrt();

    // dimensionless relative RMSE
    let relative_rmse = if ss_reference > 0.0 {
        (ss_complex / ss_reference).sqrt()
    } else {
        f64::NAN
    };

    let relative_rmse_percent = relative_rmse * 100.0;

    Ok(FitMetrics {
        rmse_real,
        rmse_imag,
        rmse_complex: rmse_magnitude,
        rmse_magnitude,
        relative_rmse,
        relative_rmse_percent,
    })
}
