use crate::drt::DrtResult;
use crate::types::EisData;
use anyhow::{Context, Result, bail};
use num_complex::Complex;
use serde::Serialize;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Serialize)]
pub struct MatlabComparison {
    pub gamma_points_compared: Option<usize>,
    pub gamma_rmse: Option<f64>,
    pub gamma_relative_rmse_percent: Option<f64>,
    pub impedance_points_compared: Option<usize>,
    pub impedance_rmse_magnitude: Option<f64>,
    pub impedance_relative_rmse_percent: Option<f64>,
    pub note: String,
}

#[derive(Debug, Clone)]
struct MatlabGammaPoint {
    tau: f64,
    gamma: f64,
}

#[derive(Debug, Clone)]
struct MatlabRegressionPoint {
    freq: f64,
    z_fit: Complex<f64>,
}

pub fn compare_with_matlab_outputs(
    data: &EisData,
    result: &DrtResult,
    matlab_drt_path: Option<&Path>,
    matlab_regression_path: Option<&Path>,
) -> Result<MatlabComparison> {
    let gamma_comparison = if let Some(path) = matlab_drt_path {
        Some(compare_gamma(result, &read_matlab_drt(path)?)?)
    } else {
        None
    };

    let impedance_comparison = if let Some(path) = matlab_regression_path {
        Some(compare_impedance(
            data,
            result,
            &read_matlab_regression(path)?,
        )?)
    } else {
        None
    };

    Ok(MatlabComparison {
        gamma_points_compared: gamma_comparison.map(|value| value.0),
        gamma_rmse: gamma_comparison.map(|value| value.1),
        gamma_relative_rmse_percent: gamma_comparison.map(|value| value.2),
        impedance_points_compared: impedance_comparison.map(|value| value.0),
        impedance_rmse_magnitude: impedance_comparison.map(|value| value.1),
        impedance_relative_rmse_percent: impedance_comparison.map(|value| value.2),
        note: "comparison against exported MATLAB DRTtools files; differences are expected unless discretization, lambda, RBF/QP settings, and preprocessing match".to_string(),
    })
}

fn read_matlab_drt(path: &Path) -> Result<Vec<MatlabGammaPoint>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read MATLAB DRT file {}", path.display()))?;
    let mut points = Vec::new();
    let mut in_data = false;
    for line in content.lines() {
        let parts: Vec<&str> = line.split(',').map(|part| part.trim()).collect();
        if parts.len() < 2 {
            continue;
        }
        if parts[0].eq_ignore_ascii_case("tau") {
            in_data = true;
            continue;
        }
        if !in_data {
            continue;
        }
        let tau = parts[0].parse::<f64>().ok();
        let gamma = parts[1].parse::<f64>().ok();
        if let (Some(tau), Some(gamma)) = (tau, gamma) {
            if tau.is_finite() && tau > 0.0 && gamma.is_finite() {
                points.push(MatlabGammaPoint { tau, gamma });
            }
        }
    }
    if points.is_empty() {
        bail!("no tau/gamma rows found in {}", path.display());
    }
    Ok(points)
}

fn read_matlab_regression(path: &Path) -> Result<Vec<MatlabRegressionPoint>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read MATLAB regression file {}", path.display()))?;
    let mut points = Vec::new();
    for line in content.lines().skip(1) {
        let parts: Vec<&str> = line.split(',').map(|part| part.trim()).collect();
        if parts.len() < 3 {
            continue;
        }
        let freq = parts[0].parse::<f64>().ok();
        let re = parts[1].parse::<f64>().ok();
        let im = parts[2].parse::<f64>().ok();
        if let (Some(freq), Some(re), Some(im)) = (freq, re, im) {
            if freq.is_finite() && freq > 0.0 && re.is_finite() && im.is_finite() {
                points.push(MatlabRegressionPoint {
                    freq,
                    z_fit: Complex::new(re, im),
                });
            }
        }
    }
    if points.is_empty() {
        bail!("no regression rows found in {}", path.display());
    }
    Ok(points)
}

fn compare_gamma(result: &DrtResult, matlab: &[MatlabGammaPoint]) -> Result<(usize, f64, f64)> {
    let mut ss = 0.0;
    let mut ss_ref = 0.0;
    let mut n = 0;
    for point in matlab {
        if let Some(gamma) = interpolate_log_tau(&result.tau, &result.gamma, point.tau) {
            let residual = gamma - point.gamma;
            ss += residual * residual;
            ss_ref += point.gamma * point.gamma;
            n += 1;
        }
    }
    if n == 0 {
        bail!("no overlapping tau range for MATLAB DRT comparison");
    }
    let rmse = (ss / n as f64).sqrt();
    let rel = (ss / ss_ref.max(1.0e-300)).sqrt() * 100.0;
    Ok((n, rmse, rel))
}

fn compare_impedance(
    data: &EisData,
    result: &DrtResult,
    matlab: &[MatlabRegressionPoint],
) -> Result<(usize, f64, f64)> {
    let mut ss = 0.0;
    let mut ss_ref = 0.0;
    let mut n = 0;
    for point in matlab {
        if let Some(z) = nearest_frequency_fit(data, result, point.freq) {
            let residual = z - point.z_fit;
            ss += residual.norm_sqr();
            ss_ref += point.z_fit.norm_sqr();
            n += 1;
        }
    }
    if n == 0 {
        bail!("no overlapping frequency range for MATLAB regression comparison");
    }
    let rmse = (ss / n as f64).sqrt();
    let rel = (ss / ss_ref.max(1.0e-300)).sqrt() * 100.0;
    Ok((n, rmse, rel))
}

fn interpolate_log_tau(tau: &[f64], gamma: &[f64], target_tau: f64) -> Option<f64> {
    if tau.len() != gamma.len() || tau.is_empty() || target_tau <= 0.0 {
        return None;
    }
    let target = target_tau.ln();
    let logs: Vec<f64> = tau.iter().map(|value| value.ln()).collect();
    if target < logs[0] || target > *logs.last()? {
        return None;
    }
    for idx in 0..logs.len().saturating_sub(1) {
        if target >= logs[idx] && target <= logs[idx + 1] {
            let span = logs[idx + 1] - logs[idx];
            let t = if span.abs() <= f64::EPSILON {
                0.0
            } else {
                (target - logs[idx]) / span
            };
            return Some(gamma[idx] + t * (gamma[idx + 1] - gamma[idx]));
        }
    }
    gamma.last().copied()
}

fn nearest_frequency_fit(data: &EisData, result: &DrtResult, freq: f64) -> Option<Complex<f64>> {
    data.frequency_hz
        .iter()
        .copied()
        .zip(result.z_fit.iter().copied())
        .min_by(|a, b| {
            (a.0.ln() - freq.ln())
                .abs()
                .total_cmp(&(b.0.ln() - freq.ln()).abs())
        })
        .map(|(_, z)| z)
}
