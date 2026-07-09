use crate::ecm::{EquivalentCircuitModel, RqrModel, RqrParams};
use crate::types::{EisData, FitMetrics, calculate_fit_metrics};
use anyhow::{Context, Result, bail};
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;
use serde::Serialize;

const N_MIN: f64 = 1.0e-6;
const N_MAX: f64 = 1.0;

#[derive(Debug, Clone, Copy, Serialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum Weighting {
    None,
    Modulus,
    Proportional,
}

#[derive(Debug, Clone)]
pub struct RqrFitSettings {
    pub initial: RqrParams,
    pub weight: Weighting,
    pub max_iter: usize,
    pub tol: f64,
}

#[derive(Debug, Clone)]
pub struct RqrFitResult {
    pub params: RqrParams,
    pub z_fit: Vec<Complex<f64>>,
    pub metrics: FitMetrics,
    pub mean_weighted_chi_square: f64,
    pub reduced_chi_square: f64,
    pub iterations: usize,
    pub converged: bool,
    pub weight: Weighting,
}

#[derive(Debug, Clone, Default)]
pub struct PartialRqrInit {
    pub rs: Option<f64>,
    pub rct: Option<f64>,
    pub q: Option<f64>,
    pub n: Option<f64>,
}

pub fn fit_rqr(data: &EisData, settings: &RqrFitSettings) -> Result<RqrFitResult> {
    settings.initial.validate()?;
    if settings.max_iter == 0 {
        bail!("max_iter must be > 0");
    }
    if settings.tol <= 0.0 || !settings.tol.is_finite() {
        bail!("tol must be finite and > 0");
    }

    let mut p = transform_params(settings.initial);
    let mut mu = 1.0e-3;

    let mut residual = residual_vector(
        data,
        inverse_transform_params(&p),
        settings.weight,
    );

    let mut cost = residual.dot(&residual);
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..settings.max_iter {
        iterations = iter + 1;

        let jacobian = finite_difference_jacobian(data, &p, settings.weight);
        let jtj = jacobian.transpose() * &jacobian;
        let jtr = jacobian.transpose() * &residual;

        let damping = DMatrix::<f64>::identity(4, 4).scale(mu);

        let step = -(jtj + damping)
            .lu()
            .solve(&jtr)
            .context("LM linear system is singular")?;

        if step.norm() < settings.tol * (p.norm() + settings.tol) {
            converged = true;
            break;
        }

        let trial_p = &p + &step;
        let trial_params = inverse_transform_params(&trial_p);

        let trial_residual = residual_vector(
            data,
            trial_params,
            settings.weight,
        );

        let trial_cost = trial_residual.dot(&trial_residual);

        if trial_cost.is_finite() && trial_cost < cost {
            let rel_improvement = (cost - trial_cost) / cost.max(1.0);

            p = trial_p;
            residual = trial_residual;
            cost = trial_cost;

            mu = (mu * 0.3).max(1.0e-12);

            if rel_improvement < settings.tol {
                converged = true;
                break;
            }
        } else {
            mu = (mu * 10.0).min(1.0e12);
        }
    }

    let params = inverse_transform_params(&p);
    let model = RqrModel { params };
    let z_fit = model.impedance(&data.frequency_hz);
    let metrics = calculate_fit_metrics(data, &z_fit)?;

    let residual_len = residual.len() as f64;
    let num_params = 4.0;
    let weighted_sse = cost;
    let mean_weighted_chi_square = if residual_len > 0.0 {
        cost / residual_len
    } else {
        f64::NAN
    };

    // reduced chi-square
    // degrees of freedom = number of data points - number of parameters
    let dof = (residual_len - num_params).max(1.0);
    let reduced_chi_square = cost / dof;

    Ok(RqrFitResult {
        params,
        z_fit,
        metrics,

        // 新增，更适合作为用户看到的 chi_squared
        mean_weighted_chi_square,
        reduced_chi_square,

        iterations,
        converged,
        weight: settings.weight,
    })
}
pub fn complete_initial_params(
    data: &EisData,
    partial: PartialRqrInit,
    auto: bool,
) -> Result<RqrParams> {
    let auto_guess = if auto {
        Some(auto_init_rqr(data))
    } else {
        None
    };

    let params = RqrParams {
        rs: partial.rs.or(auto_guess.map(|p| p.rs)).ok_or_else(|| {
            anyhow::anyhow!("missing --rs; provide all parameters or use --auto-init")
        })?,
        rct: partial.rct.or(auto_guess.map(|p| p.rct)).ok_or_else(|| {
            anyhow::anyhow!("missing --rct; provide all parameters or use --auto-init")
        })?,
        q: partial.q.or(auto_guess.map(|p| p.q)).ok_or_else(|| {
            anyhow::anyhow!("missing --q; provide all parameters or use --auto-init")
        })?,
        n: partial.n.or(auto_guess.map(|p| p.n)).ok_or_else(|| {
            anyhow::anyhow!("missing --n; provide all parameters or use --auto-init")
        })?,
    };
    params.validate()?;
    Ok(params)
}

pub fn auto_init_rqr(data: &EisData) -> RqrParams {
    let high_re = data.z_real.first().copied().unwrap_or(1.0).max(1.0e-6);
    let min_re = data
        .z_real
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .min_by(|a, b| a.total_cmp(b))
        .unwrap_or(high_re)
        .max(1.0e-6);
    let low_re = data.z_real.last().copied().unwrap_or(min_re);
    let rs = min_re.min(high_re).max(1.0e-6);
    let max_mod = data
        .z_real
        .iter()
        .zip(&data.z_imag)
        .map(|(&re, &im)| (re * re + im * im).sqrt())
        .fold(0.0_f64, f64::max);
    let rct = (low_re - rs).max(max_mod.max(1.0)).max(1.0e-6);
    let n = 0.85;

    // Rough CPE seed: use the Nyquist arc top frequency and the RC relation
    // f_peak ~= 1 / (2*pi*Rct*C), generalized here as Q ~= 1/(Rct*(2*pi*f_peak)^n).
    let f_peak = data
        .frequency_hz
        .iter()
        .zip(&data.z_imag)
        .filter(|(_, im)| im.is_finite())
        .max_by(|a, b| (-a.1).total_cmp(&(-b.1)))
        .map(|(&freq, _)| freq)
        .filter(|freq| *freq > 0.0 && freq.is_finite());
    let q = f_peak
        .map(|freq| 1.0 / (rct * (2.0 * std::f64::consts::PI * freq).powf(n)))
        .filter(|value| value.is_finite() && *value > 0.0)
        .unwrap_or(1.0e-3);

    RqrParams { rs, rct, q, n }
}

fn residual_vector(data: &EisData, params: RqrParams, weight: Weighting) -> DVector<f64> {
    let model = RqrModel { params };
    let z_fit = model.impedance(&data.frequency_hz);
    let mut values = DVector::<f64>::zeros(data.len() * 2);
    for i in 0..data.len() {
        let re = data.z_real[i];
        let im = data.z_imag[i];
        let scale = residual_scale(re, im, weight);
        values[i] = (re - z_fit[i].re) * scale;
        values[i + data.len()] = (im - z_fit[i].im) * scale;
    }
    values
}

fn residual_scale(re: f64, im: f64, weight: Weighting) -> f64 {
    match weight {
        Weighting::None => 1.0,
        Weighting::Modulus | Weighting::Proportional => {
            1.0 / (re * re + im * im).sqrt().max(1.0e-12)
        }
    }
}

fn finite_difference_jacobian(data: &EisData, p: &DVector<f64>, weight: Weighting) -> DMatrix<f64> {
    let base = residual_vector(data, inverse_transform_params(p), weight);
    let mut jacobian = DMatrix::<f64>::zeros(base.len(), p.len());
    for col in 0..p.len() {
        let step = 1.0e-6 * p[col].abs().max(1.0);
        let mut shifted = p.clone();
        shifted[col] += step;
        let shifted_residual = residual_vector(data, inverse_transform_params(&shifted), weight);
        let diff = (shifted_residual - &base).scale(1.0 / step);
        jacobian.set_column(col, &diff);
    }
    jacobian
}

fn transform_params(params: RqrParams) -> DVector<f64> {
    DVector::from_vec(vec![
        params.rs.ln(),
        params.rct.ln(),
        params.q.ln(),
        logit((params.n - N_MIN) / (N_MAX - N_MIN)),
    ])
}

fn inverse_transform_params(p: &DVector<f64>) -> RqrParams {
    RqrParams {
        rs: p[0].exp().max(1.0e-300),
        rct: p[1].exp().max(1.0e-300),
        q: p[2].exp().max(1.0e-300),
        n: N_MIN + (N_MAX - N_MIN) * sigmoid(p[3]),
    }
}

fn sigmoid(value: f64) -> f64 {
    if value >= 0.0 {
        let exp_neg = (-value).exp();
        1.0 / (1.0 + exp_neg)
    } else {
        let exp_pos = value.exp();
        exp_pos / (1.0 + exp_pos)
    }
}

fn logit(value: f64) -> f64 {
    let clipped = value.clamp(1.0e-12, 1.0 - 1.0e-12);
    (clipped / (1.0 - clipped)).ln()
}
