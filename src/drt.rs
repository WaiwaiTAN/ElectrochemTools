use crate::regularization::solve_tikhonov;
use crate::types::{EisData, FitMetrics, calculate_fit_metrics};
use anyhow::{Result, bail};
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;
use serde::Serialize;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct DrtSettings {
    pub lambda: f64,
    pub tau_min: Option<f64>,
    pub tau_max: Option<f64>,
    pub n_tau: usize,
    pub regularization_order: usize,
}

#[derive(Debug, Clone)]
pub struct DrtResult {
    pub tau: Vec<f64>,
    pub gamma: Vec<f64>,
    pub r_inf: f64,
    pub z_fit: Vec<Complex<f64>>,
    pub settings_used: DrtSettingsUsed,
    pub metrics: FitMetrics,
}

#[derive(Debug, Clone, Serialize)]
pub struct DrtSettingsUsed {
    pub lambda: f64,
    pub tau_min: f64,
    pub tau_max: f64,
    pub n_tau: usize,
    pub regularization_order: usize,
}

pub fn make_log_tau_grid(tau_min: f64, tau_max: f64, n_tau: usize) -> Result<Vec<f64>> {
    if !tau_min.is_finite() || tau_min <= 0.0 {
        bail!("tau_min must be finite and > 0");
    }
    if !tau_max.is_finite() || tau_max <= tau_min {
        bail!("tau_max must be finite and > tau_min");
    }
    if n_tau < 3 {
        bail!("n_tau must be at least 3");
    }
    let ln_min = tau_min.ln();
    let ln_max = tau_max.ln();
    let step = (ln_max - ln_min) / (n_tau as f64 - 1.0);
    Ok((0..n_tau)
        .map(|idx| (ln_min + step * idx as f64).exp())
        .collect())
}

pub fn infer_tau_bounds(data: &EisData) -> Result<(f64, f64)> {
    let f_min = data
        .frequency_hz
        .iter()
        .copied()
        .filter(|f| *f > 0.0 && f.is_finite())
        .min_by(|a, b| a.total_cmp(b))
        .ok_or_else(|| anyhow::anyhow!("cannot infer tau bounds without positive frequencies"))?;
    let f_max = data
        .frequency_hz
        .iter()
        .copied()
        .filter(|f| *f > 0.0 && f.is_finite())
        .max_by(|a, b| a.total_cmp(b))
        .ok_or_else(|| anyhow::anyhow!("cannot infer tau bounds without positive frequencies"))?;
    Ok((
        1.0 / (2.0 * PI * f_max) / 10.0,
        1.0 / (2.0 * PI * f_min) * 10.0,
    ))
}

pub fn solve_drt(data: &EisData, settings: &DrtSettings) -> Result<DrtResult> {
    if settings.lambda < 0.0 || !settings.lambda.is_finite() {
        bail!("lambda must be finite and non-negative");
    }
    let (default_min, default_max) = infer_tau_bounds(data)?;
    let tau_min = settings.tau_min.unwrap_or(default_min);
    let tau_max = settings.tau_max.unwrap_or(default_max);
    let tau = make_log_tau_grid(tau_min, tau_max, settings.n_tau)?;
    let delta_ln_tau = delta_ln_tau(&tau);
    let n_points = data.len();
    let n_tau = tau.len();
    let mut a = DMatrix::<f64>::zeros(2 * n_points, n_tau + 1);
    let mut b = DVector::<f64>::zeros(2 * n_points);

    for (i, (&freq, (&z_re, &z_im))) in data
        .frequency_hz
        .iter()
        .zip(data.z_real.iter().zip(&data.z_imag))
        .enumerate()
    {
        let omega = 2.0 * PI * freq;
        a[(i, 0)] = 1.0;
        b[i] = z_re;
        b[i + n_points] = z_im;
        for (k, (&tau_k, &dln)) in tau.iter().zip(&delta_ln_tau).enumerate() {
            let wt = omega * tau_k;
            let denom = 1.0 + wt * wt;
            a[(i, k + 1)] = dln / denom;
            a[(i + n_points, k + 1)] = -dln * wt / denom;
        }
    }

    let x = solve_tikhonov(
        &a,
        &b,
        settings.lambda,
        n_tau,
        settings.regularization_order,
    )?;
    let r_inf = x[0];
    let gamma: Vec<f64> = x.iter().skip(1).copied().collect();
    let z_fit = reconstruct_impedance(&data.frequency_hz, r_inf, &tau, &gamma);
    let metrics = calculate_fit_metrics(data, &z_fit)?;

    Ok(DrtResult {
        tau,
        gamma,
        r_inf,
        z_fit,
        settings_used: DrtSettingsUsed {
            lambda: settings.lambda,
            tau_min,
            tau_max,
            n_tau,
            regularization_order: settings.regularization_order,
        },
        metrics,
    })
}

pub fn reconstruct_impedance(
    frequency_hz: &[f64],
    r_inf: f64,
    tau: &[f64],
    gamma: &[f64],
) -> Vec<Complex<f64>> {
    let delta = delta_ln_tau(tau);
    frequency_hz
        .iter()
        .map(|&freq| {
            let omega = 2.0 * PI * freq;
            tau.iter().zip(gamma).zip(&delta).fold(
                Complex::new(r_inf, 0.0),
                |acc, ((&tau_k, &gamma_k), &dln)| {
                    let wt = omega * tau_k;
                    let denom = 1.0 + wt * wt;
                    acc + Complex::new(gamma_k * dln / denom, -gamma_k * dln * wt / denom)
                },
            )
        })
        .collect()
}

pub fn delta_ln_tau(tau: &[f64]) -> Vec<f64> {
    if tau.len() < 2 {
        return vec![1.0; tau.len()];
    }
    let logs: Vec<f64> = tau.iter().map(|value| value.ln()).collect();
    (0..tau.len())
        .map(|idx| {
            if idx == 0 {
                logs[1] - logs[0]
            } else if idx + 1 == tau.len() {
                logs[idx] - logs[idx - 1]
            } else {
                0.5 * (logs[idx + 1] - logs[idx - 1])
            }
        })
        .collect()
}
