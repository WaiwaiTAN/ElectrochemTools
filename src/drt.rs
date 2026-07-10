use crate::regularization::{
    solve_tikhonov_active_set_with_penalty, solve_tikhonov_general_with_penalty,
};
use crate::types::{EisData, FitMetrics, calculate_fit_metrics};
use anyhow::{Result, bail};
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;
use serde::Serialize;
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, Serialize, clap::ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum TauGridMode {
    Logspace,
    Drttools,
}

#[derive(Debug, Clone)]
pub struct DrtSettings {
    pub lambda: f64,
    pub tau_min: Option<f64>,
    pub tau_max: Option<f64>,
    pub n_tau: usize,
    pub tau_grid: TauGridMode,
    pub fit_inductance: bool,
    pub regularization_order: usize,
    pub nonnegative: bool,
    pub credible_intervals: bool,
}

#[derive(Debug, Clone)]
pub struct DrtResult {
    pub tau: Vec<f64>,
    pub gamma: Vec<f64>,
    pub r_inf: f64,
    pub z_fit: Vec<Complex<f64>>,
    pub settings_used: DrtSettingsUsed,
    pub metrics: FitMetrics,
    pub polarization_resistance: f64,
    pub inductance: f64,
    pub peaks: Vec<DrtPeak>,
    pub credible_intervals: Option<DrtCredibleIntervals>,
    pub kk: KkConsistencyResult,
}

#[derive(Debug, Clone, Serialize)]
pub struct DrtSettingsUsed {
    pub lambda: f64,
    pub tau_min: f64,
    pub tau_max: f64,
    pub n_tau: usize,
    pub tau_grid: TauGridMode,
    pub fit_inductance: bool,
    pub regularization_order: usize,
    pub nonnegative: bool,
    pub credible_intervals: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct DrtPeak {
    pub tau: f64,
    pub log10_tau: f64,
    pub gamma: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct LambdaScanPoint {
    pub lambda: f64,
    pub relative_rmse: f64,
    pub roughness: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct DrtCredibleIntervals {
    pub inductance_std: Option<f64>,
    pub r_inf_std: f64,
    pub gamma_std: Vec<f64>,
    pub gamma_lower_95: Vec<f64>,
    pub gamma_upper_95: Vec<f64>,
    pub note: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct KkConsistencyResult {
    pub real_to_imag_relative_rmse_percent: f64,
    pub imag_to_real_relative_rmse_percent: f64,
    pub mean_score: f64,
    pub note: String,
    #[serde(skip_serializing)]
    pub rows: Vec<KkConsistencyRow>,
}

#[derive(Debug, Clone, Serialize)]
pub struct KkConsistencyRow {
    pub frequency: f64,
    pub z_real_exp: f64,
    pub z_imag_exp: f64,
    pub z_real_from_imag: f64,
    pub z_imag_from_real: f64,
    pub residual_real_from_imag: f64,
    pub residual_imag_from_real: f64,
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

pub fn make_drttools_tau_grid(data: &EisData) -> Result<Vec<f64>> {
    let mut tau: Vec<f64> = data
        .frequency_hz
        .iter()
        .copied()
        .filter(|freq| *freq > 0.0 && freq.is_finite())
        .map(|freq| 1.0 / freq)
        .collect();
    if tau.len() < 3 {
        bail!("DRTtools-style tau grid requires at least three positive frequencies");
    }
    tau.sort_by(|a, b| a.total_cmp(b));
    tau.dedup_by(|a, b| (*a - *b).abs() <= 1.0e-15 * a.abs().max(b.abs()).max(1.0));
    Ok(tau)
}

pub fn solve_drt(data: &EisData, settings: &DrtSettings) -> Result<DrtResult> {
    if settings.lambda < 0.0 || !settings.lambda.is_finite() {
        bail!("lambda must be finite and non-negative");
    }
    let (default_min, default_max) = infer_tau_bounds(data)?;
    let tau_min = settings.tau_min.unwrap_or(default_min);
    let tau_max = settings.tau_max.unwrap_or(default_max);
    let tau = match settings.tau_grid {
        TauGridMode::Logspace => make_log_tau_grid(tau_min, tau_max, settings.n_tau)?,
        TauGridMode::Drttools => make_drttools_tau_grid(data)?,
    };
    let n_tau = tau.len();
    let n_unpenalized = if settings.fit_inductance { 2 } else { 1 };
    let gamma_offset = n_unpenalized;
    let (a, b) = build_drt_system(data, &tau, settings.fit_inductance);
    let penalty =
        drttools_piecewise_linear_penalty(&tau, settings.regularization_order, n_unpenalized)?;

    let x = if settings.nonnegative {
        let lower_bound_zero = vec![true; n_unpenalized + n_tau];
        solve_tikhonov_active_set_with_penalty(
            &a,
            &b,
            settings.lambda,
            &penalty,
            &lower_bound_zero,
            1_000,
            1.0e-9,
        )?
    } else {
        solve_tikhonov_general_with_penalty(&a, &b, settings.lambda, &penalty)?
    };
    let (inductance, r_inf) = if settings.fit_inductance {
        (x[0], x[1])
    } else {
        (0.0, x[0])
    };
    let gamma: Vec<f64> = x.iter().skip(gamma_offset).copied().collect();
    let z_fit =
        reconstruct_impedance_with_inductance(&data.frequency_hz, r_inf, inductance, &tau, &gamma);
    let metrics = calculate_fit_metrics(data, &z_fit)?;
    let polarization_resistance = polarization_resistance(&tau, &gamma);
    let peaks = detect_drt_peaks(&tau, &gamma, 0.01);
    let used_tau_min = tau.first().copied().unwrap_or(tau_min);
    let used_tau_max = tau.last().copied().unwrap_or(tau_max);
    let credible_intervals = if settings.credible_intervals {
        estimate_drt_credible_intervals(
            &a,
            &b,
            &x,
            n_tau,
            settings.lambda,
            &penalty,
            n_unpenalized,
            settings.fit_inductance,
        )
    } else {
        None
    };
    let kk = estimate_kk_consistency(data, &tau, settings.lambda, settings.regularization_order)?;

    Ok(DrtResult {
        tau,
        gamma,
        r_inf,
        z_fit,
        settings_used: DrtSettingsUsed {
            lambda: settings.lambda,
            tau_min: used_tau_min,
            tau_max: used_tau_max,
            n_tau,
            tau_grid: settings.tau_grid,
            fit_inductance: settings.fit_inductance,
            regularization_order: settings.regularization_order,
            nonnegative: settings.nonnegative,
            credible_intervals: settings.credible_intervals,
        },
        metrics,
        polarization_resistance,
        inductance,
        peaks,
        credible_intervals,
        kk,
    })
}

pub fn scan_lambda(
    data: &EisData,
    settings: &DrtSettings,
    lambda_min: f64,
    lambda_max: f64,
    n_lambda: usize,
) -> Result<(f64, Vec<LambdaScanPoint>)> {
    if lambda_min <= 0.0 || lambda_max <= lambda_min || n_lambda < 2 {
        bail!("auto-lambda scan requires 0 < lambda_min < lambda_max and n_lambda >= 2");
    }
    let log_min = lambda_min.log10();
    let log_max = lambda_max.log10();
    let mut scan = Vec::new();
    for idx in 0..n_lambda {
        let t = idx as f64 / (n_lambda as f64 - 1.0);
        let lambda = 10_f64.powf(log_min + t * (log_max - log_min));
        let mut trial = settings.clone();
        trial.lambda = lambda;
        trial.credible_intervals = false;
        let result = solve_drt(data, &trial)?;
        let roughness = gamma_roughness(&result.gamma);
        let score = result.metrics.relative_rmse + 1.0e-3 * roughness.ln_1p();
        scan.push(LambdaScanPoint {
            lambda,
            relative_rmse: result.metrics.relative_rmse,
            roughness,
            score,
        });
    }
    let best_lambda = scan
        .iter()
        .min_by(|a, b| a.score.total_cmp(&b.score))
        .map(|point| point.lambda)
        .ok_or_else(|| anyhow::anyhow!("lambda scan produced no candidates"))?;
    Ok((best_lambda, scan))
}

pub fn detect_drt_peaks(tau: &[f64], gamma: &[f64], min_relative_height: f64) -> Vec<DrtPeak> {
    if tau.len() != gamma.len() || tau.len() < 3 {
        return Vec::new();
    }
    let max_gamma = gamma
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    if !max_gamma.is_finite() || max_gamma <= 0.0 {
        return Vec::new();
    }
    let threshold = max_gamma * min_relative_height.max(0.0);
    let mut peaks = Vec::new();
    if gamma[0] >= threshold && gamma[0] >= gamma[1] {
        peaks.push(DrtPeak {
            tau: tau[0],
            log10_tau: tau[0].log10(),
            gamma: gamma[0],
        });
    }
    for idx in 1..(gamma.len() - 1) {
        if gamma[idx] >= threshold && gamma[idx] >= gamma[idx - 1] && gamma[idx] >= gamma[idx + 1] {
            peaks.push(DrtPeak {
                tau: tau[idx],
                log10_tau: tau[idx].log10(),
                gamma: gamma[idx],
            });
        }
    }
    let last = gamma.len() - 1;
    if gamma[last] >= threshold && gamma[last] >= gamma[last - 1] {
        peaks.push(DrtPeak {
            tau: tau[last],
            log10_tau: tau[last].log10(),
            gamma: gamma[last],
        });
    }
    peaks.sort_by(|a, b| b.gamma.total_cmp(&a.gamma));
    peaks
}

fn build_drt_system(
    data: &EisData,
    tau: &[f64],
    fit_inductance: bool,
) -> (DMatrix<f64>, DVector<f64>) {
    let delta_ln_tau = delta_ln_tau(tau);
    let n_points = data.len();
    let n_tau = tau.len();
    let n_unpenalized = if fit_inductance { 2 } else { 1 };
    let gamma_offset = n_unpenalized;
    let r_col = if fit_inductance { 1 } else { 0 };
    let mut a = DMatrix::<f64>::zeros(2 * n_points, n_tau + n_unpenalized);
    let mut b = DVector::<f64>::zeros(2 * n_points);

    for (i, (&freq, (&z_re, &z_im))) in data
        .frequency_hz
        .iter()
        .zip(data.z_real.iter().zip(&data.z_imag))
        .enumerate()
    {
        let omega = 2.0 * PI * freq;
        a[(i, r_col)] = 1.0;
        if fit_inductance {
            a[(i + n_points, 0)] = omega;
        }
        b[i] = z_re;
        b[i + n_points] = z_im;
        for (k, (&tau_k, &dln)) in tau.iter().zip(&delta_ln_tau).enumerate() {
            let wt = omega * tau_k;
            let denom = 1.0 + wt * wt;
            a[(i, k + gamma_offset)] = dln / denom;
            a[(i + n_points, k + gamma_offset)] = -dln * wt / denom;
        }
    }
    (a, b)
}

fn build_real_system(data: &EisData, tau: &[f64]) -> (DMatrix<f64>, DVector<f64>) {
    let delta_ln_tau = delta_ln_tau(tau);
    let mut a = DMatrix::<f64>::zeros(data.len(), tau.len() + 1);
    let mut b = DVector::<f64>::zeros(data.len());
    for (i, (&freq, &z_re)) in data.frequency_hz.iter().zip(&data.z_real).enumerate() {
        let omega = 2.0 * PI * freq;
        a[(i, 0)] = 1.0;
        b[i] = z_re;
        for (k, (&tau_k, &dln)) in tau.iter().zip(&delta_ln_tau).enumerate() {
            let wt = omega * tau_k;
            a[(i, k + 1)] = dln / (1.0 + wt * wt);
        }
    }
    (a, b)
}

fn build_imag_system(data: &EisData, tau: &[f64]) -> (DMatrix<f64>, DVector<f64>) {
    let delta_ln_tau = delta_ln_tau(tau);
    let mut a = DMatrix::<f64>::zeros(data.len(), tau.len() + 1);
    let mut b = DVector::<f64>::zeros(data.len());
    for (i, (&freq, &z_im)) in data.frequency_hz.iter().zip(&data.z_imag).enumerate() {
        let omega = 2.0 * PI * freq;
        b[i] = z_im;
        for (k, (&tau_k, &dln)) in tau.iter().zip(&delta_ln_tau).enumerate() {
            let wt = omega * tau_k;
            a[(i, k + 1)] = -dln * wt / (1.0 + wt * wt);
        }
    }
    (a, b)
}

#[allow(clippy::too_many_arguments)]
fn estimate_drt_credible_intervals(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    x: &DVector<f64>,
    n_gamma: usize,
    lambda: f64,
    penalty: &DMatrix<f64>,
    n_unpenalized: usize,
    fit_inductance: bool,
) -> Option<DrtCredibleIntervals> {
    let residual = a * x - b;
    let dof = (b.len() as f64 - x.len() as f64).max(1.0);
    let sigma2 = residual.dot(&residual) / dof;
    let precision = a.transpose() * a + penalty.scale(lambda);
    let covariance = precision.try_inverse()?.scale(sigma2);
    let mut gamma_std = Vec::with_capacity(n_gamma);
    let mut gamma_lower_95 = Vec::with_capacity(n_gamma);
    let mut gamma_upper_95 = Vec::with_capacity(n_gamma);
    for idx in 0..n_gamma {
        let param_idx = idx + n_unpenalized;
        let std = covariance[(param_idx, param_idx)].max(0.0).sqrt();
        let center = x[param_idx];
        gamma_std.push(std);
        gamma_lower_95.push(center - 1.96 * std);
        gamma_upper_95.push(center + 1.96 * std);
    }
    let r_inf_idx = if fit_inductance { 1 } else { 0 };
    Some(DrtCredibleIntervals {
        inductance_std: fit_inductance.then(|| covariance[(0, 0)].max(0.0).sqrt()),
        r_inf_std: covariance[(r_inf_idx, r_inf_idx)].max(0.0).sqrt(),
        gamma_std,
        gamma_lower_95,
        gamma_upper_95,
        note:
            "linear-Gaussian posterior approximation around the Tikhonov solution; not HMC sampling"
                .to_string(),
    })
}

pub fn estimate_kk_consistency(
    data: &EisData,
    tau: &[f64],
    lambda: f64,
    order: usize,
) -> Result<KkConsistencyResult> {
    let (a_re, b_re) = build_real_system(data, tau);
    let (a_im, b_im) = build_imag_system(data, tau);
    let penalty = drttools_piecewise_linear_penalty(tau, order, 1)?;
    let x_re = solve_tikhonov_general_with_penalty(&a_re, &b_re, lambda, &penalty)?;
    let x_im = solve_tikhonov_general_with_penalty(&a_im, &b_im, lambda, &penalty)?;
    let gamma_re: Vec<f64> = x_re.iter().skip(1).copied().collect();
    let gamma_im: Vec<f64> = x_im.iter().skip(1).copied().collect();

    let z_from_re = reconstruct_impedance(&data.frequency_hz, x_re[0], tau, &gamma_re);
    let z_real_no_offset = reconstruct_impedance(&data.frequency_hz, 0.0, tau, &gamma_im);
    let r_inf_from_imag = data
        .z_real
        .iter()
        .zip(&z_real_no_offset)
        .map(|(&exp, fit)| exp - fit.re)
        .sum::<f64>()
        / data.len() as f64;
    let z_from_im = reconstruct_impedance(&data.frequency_hz, r_inf_from_imag, tau, &gamma_im);

    let mut ss_im = 0.0;
    let mut ref_im = 0.0;
    let mut ss_re = 0.0;
    let mut ref_re = 0.0;
    let mut rows = Vec::with_capacity(data.len());
    for i in 0..data.len() {
        let residual_im = data.z_imag[i] - z_from_re[i].im;
        let residual_re = data.z_real[i] - z_from_im[i].re;
        ss_im += residual_im * residual_im;
        ref_im += data.z_imag[i] * data.z_imag[i];
        ss_re += residual_re * residual_re;
        ref_re += data.z_real[i] * data.z_real[i];
        rows.push(KkConsistencyRow {
            frequency: data.frequency_hz[i],
            z_real_exp: data.z_real[i],
            z_imag_exp: data.z_imag[i],
            z_real_from_imag: z_from_im[i].re,
            z_imag_from_real: z_from_re[i].im,
            residual_real_from_imag: residual_re,
            residual_imag_from_real: residual_im,
        });
    }

    let real_to_imag = (ss_im / ref_im.max(1.0e-300)).sqrt();
    let imag_to_real = (ss_re / ref_re.max(1.0e-300)).sqrt();
    let mean_score =
        0.5 * ((1.0 - real_to_imag).clamp(0.0, 1.0) + (1.0 - imag_to_real).clamp(0.0, 1.0));

    Ok(KkConsistencyResult {
        real_to_imag_relative_rmse_percent: real_to_imag * 100.0,
        imag_to_real_relative_rmse_percent: imag_to_real * 100.0,
        mean_score,
        note: "linear DRT-based Hilbert/Kramers-Kronig consistency proxy; not the full Bayesian DRTtools score"
            .to_string(),
        rows,
    })
}

fn gamma_roughness(gamma: &[f64]) -> f64 {
    gamma
        .windows(2)
        .map(|pair| {
            let delta = pair[1] - pair[0];
            delta * delta
        })
        .sum::<f64>()
        .sqrt()
}

fn polarization_resistance(tau: &[f64], gamma: &[f64]) -> f64 {
    let delta = delta_ln_tau(tau);
    gamma
        .iter()
        .zip(delta)
        .map(|(&gamma, dln)| gamma * dln)
        .sum()
}

fn drttools_piecewise_linear_penalty(
    tau: &[f64],
    order: usize,
    n_unpenalized: usize,
) -> Result<DMatrix<f64>> {
    let n_gamma = tau.len();
    let n_params = n_gamma + n_unpenalized;
    let rows = match order {
        0 => n_gamma,
        1 => {
            if n_gamma < 2 {
                bail!("first-order regularization requires at least two gamma values");
            }
            n_gamma - 1
        }
        2 => {
            if n_gamma < 3 {
                bail!("second-order regularization requires at least three gamma values");
            }
            n_gamma - 2
        }
        _ => bail!("regularization order must be 0, 1, or 2"),
    };
    let mut l = DMatrix::<f64>::zeros(rows, n_params);
    match order {
        0 => {
            for row in 0..n_gamma {
                l[(row, row + n_unpenalized)] = 1.0;
            }
        }
        1 => {
            for row in 0..(n_gamma - 1) {
                let delta = (tau[row + 1] / tau[row]).ln();
                if delta <= 0.0 || !delta.is_finite() {
                    bail!("tau grid must be strictly increasing for DRT regularization");
                }
                l[(row, row + n_unpenalized)] = -1.0 / delta;
                l[(row, row + n_unpenalized + 1)] = 1.0 / delta;
            }
        }
        2 => {
            for row in 0..(n_gamma - 2) {
                let delta = (tau[row + 1] / tau[row]).ln();
                if delta <= 0.0 || !delta.is_finite() {
                    bail!("tau grid must be strictly increasing for DRT regularization");
                }
                let scale = if row == 0 || row + 1 == n_gamma - 2 {
                    2.0 / (delta * delta)
                } else {
                    1.0 / (delta * delta)
                };
                l[(row, row + n_unpenalized)] = scale;
                l[(row, row + n_unpenalized + 1)] = -2.0 * scale;
                l[(row, row + n_unpenalized + 2)] = scale;
            }
        }
        _ => unreachable!(),
    }
    Ok(l.transpose() * l)
}

pub fn reconstruct_impedance(
    frequency_hz: &[f64],
    r_inf: f64,
    tau: &[f64],
    gamma: &[f64],
) -> Vec<Complex<f64>> {
    reconstruct_impedance_with_inductance(frequency_hz, r_inf, 0.0, tau, gamma)
}

pub fn reconstruct_impedance_with_inductance(
    frequency_hz: &[f64],
    r_inf: f64,
    inductance: f64,
    tau: &[f64],
    gamma: &[f64],
) -> Vec<Complex<f64>> {
    let delta = delta_ln_tau(tau);
    frequency_hz
        .iter()
        .map(|&freq| {
            let omega = 2.0 * PI * freq;
            tau.iter().zip(gamma).zip(&delta).fold(
                Complex::new(r_inf, omega * inductance),
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
                0.5 * (logs[1] - logs[0])
            } else if idx + 1 == tau.len() {
                0.5 * (logs[idx] - logs[idx - 1])
            } else {
                0.5 * (logs[idx + 1] - logs[idx - 1])
            }
        })
        .collect()
}
