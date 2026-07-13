use crate::ecm::{
    BranchKind, EcmModel, EcmModelSpec, EcmParams, EquivalentCircuitModel, ParallelBranchParams,
    ReactiveElementParams, RqrParams,
};
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
pub struct EcmFitSettings {
    pub model: EcmModelSpec,
    pub initial: EcmParams,
    pub weight: Weighting,
    pub max_iter: usize,
    pub tol: f64,
}

#[derive(Debug, Clone)]
pub struct EcmFitResult {
    pub model: EcmModelSpec,
    pub params: EcmParams,
    pub z_fit: Vec<Complex<f64>>,
    pub metrics: FitMetrics,
    pub weighted_sse: f64,
    pub mean_weighted_chi_square: f64,
    pub reduced_chi_square: f64,
    pub parameter_std_errors: Option<Vec<f64>>,
    pub parameter_rel_std_error_percent: Option<Vec<f64>>,
    pub parameter_correlation_matrix: Option<Vec<Vec<f64>>>,
    pub iterations: usize,
    pub converged: bool,
    pub weight: Weighting,
}

#[derive(Debug, Clone, Default)]
pub struct PartialEcmInit {
    pub rs: Option<f64>,
    pub r1: Option<f64>,
    pub c1: Option<f64>,
    pub q1: Option<f64>,
    pub n1: Option<f64>,
    pub r2: Option<f64>,
    pub c2: Option<f64>,
    pub q2: Option<f64>,
    pub n2: Option<f64>,
    pub warburg_sigma: Option<f64>,
}

pub fn fit_ecm(data: &EisData, settings: &EcmFitSettings) -> Result<EcmFitResult> {
    settings.initial.validate_for(&settings.model)?;
    if settings.max_iter == 0 {
        bail!("max_iter must be > 0");
    }
    if settings.tol <= 0.0 || !settings.tol.is_finite() {
        bail!("tol must be finite and > 0");
    }

    let mut transformed = transform_params(&settings.model, &settings.initial);
    let mut damping_factor = 1.0e-3;
    let mut residual = residual_vector(data, &settings.model, &transformed, settings.weight);
    let mut cost = residual.dot(&residual);
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..settings.max_iter {
        iterations = iter + 1;
        let jacobian =
            finite_difference_jacobian(data, &settings.model, &transformed, settings.weight);
        let jtj = jacobian.transpose() * &jacobian;
        let jtr = jacobian.transpose() * &residual;
        let damping =
            DMatrix::<f64>::identity(transformed.len(), transformed.len()).scale(damping_factor);
        let step = -(jtj + damping)
            .lu()
            .solve(&jtr)
            .context("LM linear system is singular")?;

        if step.norm() < settings.tol * (transformed.norm() + settings.tol) {
            converged = true;
            break;
        }

        let trial = &transformed + &step;
        let trial_residual = residual_vector(data, &settings.model, &trial, settings.weight);
        let trial_cost = trial_residual.dot(&trial_residual);
        if trial_cost.is_finite() && trial_cost < cost {
            let relative_improvement = (cost - trial_cost) / cost.max(1.0);
            transformed = trial;
            residual = trial_residual;
            cost = trial_cost;
            damping_factor = (damping_factor * 0.3).max(1.0e-12);
            if relative_improvement < settings.tol {
                converged = true;
                break;
            }
        } else {
            damping_factor = (damping_factor * 10.0).min(1.0e12);
        }
    }

    let params = inverse_transform_params(&settings.model, &transformed);
    let fitted_model = EcmModel {
        spec: settings.model.clone(),
        params: params.clone(),
    };
    let z_fit = fitted_model.impedance(&data.frequency_hz);
    let metrics = calculate_fit_metrics(data, &z_fit)?;
    let residual_len = residual.len() as f64;
    let weighted_sse = cost;
    let mean_weighted_chi_square = cost / residual_len;
    let degrees_of_freedom = (residual_len - transformed.len() as f64).max(1.0);
    let reduced_chi_square = cost / degrees_of_freedom;
    let uncertainty = estimate_uncertainty(
        data,
        &settings.model,
        &transformed,
        settings.weight,
        reduced_chi_square,
    );

    Ok(EcmFitResult {
        model: settings.model.clone(),
        params,
        z_fit,
        metrics,
        weighted_sse,
        mean_weighted_chi_square,
        reduced_chi_square,
        parameter_std_errors: uncertainty.as_ref().map(|value| value.std_errors.clone()),
        parameter_rel_std_error_percent: uncertainty
            .as_ref()
            .map(|value| value.rel_std_error_percent.clone()),
        parameter_correlation_matrix: uncertainty.map(|value| value.correlation_matrix),
        iterations,
        converged,
        weight: settings.weight,
    })
}

pub fn complete_initial_ecm(
    data: &EisData,
    model: &EcmModelSpec,
    partial: PartialEcmInit,
    auto: bool,
) -> Result<EcmParams> {
    validate_relevant_initial_flags(model, &partial)?;
    let automatic = auto.then(|| auto_init_ecm(data, model));
    let rs = choose_initial("--rs", partial.rs, automatic.as_ref().map(|p| p.rs))?;
    let mut branches = Vec::with_capacity(model.branches.len());
    for (index, kind) in model.branches.iter().enumerate() {
        let number = index + 1;
        let auto_branch = automatic.as_ref().map(|p| p.branches[index]);
        let (r, c, q, n) = if index == 0 {
            (partial.r1, partial.c1, partial.q1, partial.n1)
        } else {
            (partial.r2, partial.c2, partial.q2, partial.n2)
        };
        let resistance = choose_initial(
            &format!("--r{number}"),
            r,
            auto_branch.map(|branch| branch.r),
        )?;
        let reactive = match kind {
            BranchKind::Rc => ReactiveElementParams::Capacitor {
                c: choose_initial(
                    &format!("--c{number}"),
                    c,
                    auto_branch.and_then(|branch| match branch.reactive {
                        ReactiveElementParams::Capacitor { c } => Some(c),
                        _ => None,
                    }),
                )?,
            },
            BranchKind::Rq => ReactiveElementParams::Cpe {
                q: choose_initial(
                    &format!("--q{number}"),
                    q,
                    auto_branch.and_then(|branch| match branch.reactive {
                        ReactiveElementParams::Cpe { q, .. } => Some(q),
                        _ => None,
                    }),
                )?,
                n: choose_initial(
                    &format!("--n{number}"),
                    n,
                    auto_branch.and_then(|branch| match branch.reactive {
                        ReactiveElementParams::Cpe { n, .. } => Some(n),
                        _ => None,
                    }),
                )?,
            },
        };
        branches.push(ParallelBranchParams {
            r: resistance,
            reactive,
        });
    }
    let warburg_sigma = if model.warburg {
        Some(choose_initial(
            "--warburg-sigma",
            partial.warburg_sigma,
            automatic.as_ref().and_then(|p| p.warburg_sigma),
        )?)
    } else {
        None
    };
    let params = EcmParams {
        rs,
        branches,
        warburg_sigma,
    };
    params.validate_for(model)?;
    Ok(params)
}

pub fn auto_init_ecm(data: &EisData, model: &EcmModelSpec) -> EcmParams {
    let high_frequency_index = data
        .frequency_hz
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index)
        .unwrap_or(0);
    let low_frequency_index = data
        .frequency_hz
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index)
        .unwrap_or(0);
    let minimum_real = data
        .z_real
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .min_by(f64::total_cmp)
        .unwrap_or(1.0)
        .max(1.0e-6);
    let rs = data.z_real[high_frequency_index]
        .min(minimum_real)
        .max(1.0e-6);
    let real_span = (data.z_real[low_frequency_index] - rs).max(1.0e-6);
    let resistance_per_branch = (real_span / model.branches.len() as f64).max(1.0e-6);
    let positive_frequencies = data
        .frequency_hz
        .iter()
        .copied()
        .filter(|f| f.is_finite() && *f > 0.0)
        .collect::<Vec<_>>();
    let f_min = positive_frequencies
        .iter()
        .copied()
        .min_by(f64::total_cmp)
        .unwrap_or(1.0e-3);
    let f_max = positive_frequencies
        .iter()
        .copied()
        .max_by(f64::total_cmp)
        .unwrap_or(1.0e3);
    let log_span = (f_max / f_min).max(1.0).ln();
    let branches = model
        .branches
        .iter()
        .enumerate()
        .map(|(index, kind)| {
            let fraction = (index + 1) as f64 / (model.branches.len() + 1) as f64;
            let characteristic_frequency = f_max * (-fraction * log_span).exp();
            let omega = 2.0 * std::f64::consts::PI * characteristic_frequency;
            let reactive = match kind {
                BranchKind::Rc => ReactiveElementParams::Capacitor {
                    c: 1.0 / (resistance_per_branch * omega),
                },
                BranchKind::Rq => {
                    let n = 0.85;
                    ReactiveElementParams::Cpe {
                        q: 1.0 / (resistance_per_branch * omega.powf(n)),
                        n,
                    }
                }
            };
            ParallelBranchParams {
                r: resistance_per_branch,
                reactive,
            }
        })
        .collect();
    let warburg_sigma = model.warburg.then(|| {
        let low_imag = (-data.z_imag[low_frequency_index]).max(0.1 * real_span);
        (low_imag * (2.0 * std::f64::consts::PI * f_min).sqrt()).max(1.0e-12)
    });
    EcmParams {
        rs,
        branches,
        warburg_sigma,
    }
}

fn validate_relevant_initial_flags(model: &EcmModelSpec, partial: &PartialEcmInit) -> Result<()> {
    if model.branches.len() == 1
        && [partial.r2, partial.c2, partial.q2, partial.n2]
            .iter()
            .any(Option::is_some)
    {
        bail!("second-branch parameters were supplied for a one-branch model");
    }
    for (index, kind) in model.branches.iter().enumerate() {
        let (c, q, n) = if index == 0 {
            (partial.c1, partial.q1, partial.n1)
        } else {
            (partial.c2, partial.q2, partial.n2)
        };
        match kind {
            BranchKind::Rc if q.is_some() || n.is_some() => {
                bail!("Q/n parameters were supplied for RC branch {}", index + 1)
            }
            BranchKind::Rq if c.is_some() => {
                bail!("C parameters were supplied for RQ branch {}", index + 1)
            }
            _ => {}
        }
    }
    if !model.warburg && partial.warburg_sigma.is_some() {
        bail!("--warburg-sigma requires a model name ending in _W");
    }
    Ok(())
}

fn choose_initial(label: &str, supplied: Option<f64>, automatic: Option<f64>) -> Result<f64> {
    supplied.or(automatic).ok_or_else(|| {
        anyhow::anyhow!("missing {label}; provide all parameters or use --auto-init")
    })
}

fn residual_vector(
    data: &EisData,
    model: &EcmModelSpec,
    transformed: &DVector<f64>,
    weight: Weighting,
) -> DVector<f64> {
    let fitted_model = EcmModel {
        spec: model.clone(),
        params: inverse_transform_params(model, transformed),
    };
    let z_fit = fitted_model.impedance(&data.frequency_hz);
    let mut values = DVector::<f64>::zeros(data.len() * 2);
    for index in 0..data.len() {
        let re = data.z_real[index];
        let im = data.z_imag[index];
        let scale = residual_scale(re, im, weight);
        values[index] = (re - z_fit[index].re) * scale;
        values[index + data.len()] = (im - z_fit[index].im) * scale;
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

fn finite_difference_jacobian(
    data: &EisData,
    model: &EcmModelSpec,
    transformed: &DVector<f64>,
    weight: Weighting,
) -> DMatrix<f64> {
    let base = residual_vector(data, model, transformed, weight);
    let mut jacobian = DMatrix::<f64>::zeros(base.len(), transformed.len());
    for column in 0..transformed.len() {
        let step = 1.0e-6 * transformed[column].abs().max(1.0);
        let mut shifted = transformed.clone();
        shifted[column] += step;
        let difference = (residual_vector(data, model, &shifted, weight) - &base).scale(1.0 / step);
        jacobian.set_column(column, &difference);
    }
    jacobian
}

#[derive(Debug, Clone)]
struct ParameterUncertainty {
    std_errors: Vec<f64>,
    rel_std_error_percent: Vec<f64>,
    correlation_matrix: Vec<Vec<f64>>,
}

fn estimate_uncertainty(
    data: &EisData,
    model: &EcmModelSpec,
    transformed: &DVector<f64>,
    weight: Weighting,
    reduced_chi_square: f64,
) -> Option<ParameterUncertainty> {
    let jacobian = finite_difference_jacobian(data, model, transformed, weight);
    let transformed_covariance = (jacobian.transpose() * jacobian)
        .try_inverse()?
        .scale(reduced_chi_square);
    let params = inverse_transform_params(model, transformed);
    let values = params.values();
    let derivatives = parameter_derivatives(model, &params);
    let transform = DMatrix::<f64>::from_diagonal(&DVector::from_vec(derivatives));
    let covariance = &transform * transformed_covariance * transform.transpose();
    if (0..values.len()).any(|index| {
        let variance = covariance[(index, index)];
        !variance.is_finite() || variance <= 0.0
    }) {
        return None;
    }
    let std_errors = (0..values.len())
        .map(|index| covariance[(index, index)].sqrt())
        .collect::<Vec<_>>();
    let mut correlation_matrix = vec![vec![f64::NAN; values.len()]; values.len()];
    for row in 0..values.len() {
        for column in 0..values.len() {
            let denominator = std_errors[row] * std_errors[column];
            correlation_matrix[row][column] = if denominator > 0.0 {
                (covariance[(row, column)] / denominator).clamp(-1.0, 1.0)
            } else if row == column {
                1.0
            } else {
                f64::NAN
            };
        }
    }
    let rel_std_error_percent = std_errors
        .iter()
        .zip(&values)
        .map(|(std, value)| 100.0 * std / value.abs().max(1.0e-300))
        .collect();
    Some(ParameterUncertainty {
        std_errors,
        rel_std_error_percent,
        correlation_matrix,
    })
}

fn transform_params(model: &EcmModelSpec, params: &EcmParams) -> DVector<f64> {
    let mut values = Vec::with_capacity(model.parameter_labels().len());
    values.push(params.rs.ln());
    for branch in &params.branches {
        values.push(branch.r.ln());
        match branch.reactive {
            ReactiveElementParams::Capacitor { c } => values.push(c.ln()),
            ReactiveElementParams::Cpe { q, n } => {
                values.push(q.ln());
                values.push(logit((n - N_MIN) / (N_MAX - N_MIN)));
            }
        }
    }
    if let Some(sigma) = params.warburg_sigma {
        values.push(sigma.ln());
    }
    DVector::from_vec(values)
}

fn inverse_transform_params(model: &EcmModelSpec, transformed: &DVector<f64>) -> EcmParams {
    let mut index = 0;
    let rs = transformed[index].exp().max(1.0e-300);
    index += 1;
    let mut branches = Vec::with_capacity(model.branches.len());
    for kind in &model.branches {
        let r = transformed[index].exp().max(1.0e-300);
        index += 1;
        let reactive = match kind {
            BranchKind::Rc => {
                let c = transformed[index].exp().max(1.0e-300);
                index += 1;
                ReactiveElementParams::Capacitor { c }
            }
            BranchKind::Rq => {
                let q = transformed[index].exp().max(1.0e-300);
                let n = N_MIN + (N_MAX - N_MIN) * sigmoid(transformed[index + 1]);
                index += 2;
                ReactiveElementParams::Cpe { q, n }
            }
        };
        branches.push(ParallelBranchParams { r, reactive });
    }
    let warburg_sigma = model
        .warburg
        .then(|| transformed[index].exp().max(1.0e-300));
    EcmParams {
        rs,
        branches,
        warburg_sigma,
    }
}

fn parameter_derivatives(model: &EcmModelSpec, params: &EcmParams) -> Vec<f64> {
    let mut derivatives = vec![params.rs];
    for branch in &params.branches {
        derivatives.push(branch.r);
        match branch.reactive {
            ReactiveElementParams::Capacitor { c } => derivatives.push(c),
            ReactiveElementParams::Cpe { q, n } => {
                derivatives.push(q);
                let scaled = (n - N_MIN) / (N_MAX - N_MIN);
                derivatives.push((N_MAX - N_MIN) * scaled * (1.0 - scaled));
            }
        }
    }
    if model.warburg {
        derivatives.push(params.warburg_sigma.expect("Warburg parameter must exist"));
    }
    derivatives
}

fn sigmoid(value: f64) -> f64 {
    if value >= 0.0 {
        let exp_negative = (-value).exp();
        1.0 / (1.0 + exp_negative)
    } else {
        let exp_positive = value.exp();
        exp_positive / (1.0 + exp_positive)
    }
}

fn logit(value: f64) -> f64 {
    let clipped = value.clamp(1.0e-12, 1.0 - 1.0e-12);
    (clipped / (1.0 - clipped)).ln()
}

// Backward-compatible R(QR) fitting API.
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
    pub weighted_sse: f64,
    pub mean_weighted_chi_square: f64,
    pub reduced_chi_square: f64,
    pub parameter_std_errors: Option<RqrParams>,
    pub parameter_rel_std_error_percent: Option<RqrParams>,
    pub parameter_correlation_matrix: Option<Vec<Vec<f64>>>,
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
    let model: EcmModelSpec = "R_QR".parse()?;
    let generic = fit_ecm(
        data,
        &EcmFitSettings {
            model,
            initial: rqr_to_ecm(settings.initial),
            weight: settings.weight,
            max_iter: settings.max_iter,
            tol: settings.tol,
        },
    )?;
    Ok(RqrFitResult {
        params: ecm_to_rqr(&generic.params),
        z_fit: generic.z_fit,
        metrics: generic.metrics,
        weighted_sse: generic.weighted_sse,
        mean_weighted_chi_square: generic.mean_weighted_chi_square,
        reduced_chi_square: generic.reduced_chi_square,
        parameter_std_errors: generic.parameter_std_errors.as_deref().map(vector_to_rqr),
        parameter_rel_std_error_percent: generic
            .parameter_rel_std_error_percent
            .as_deref()
            .map(vector_to_rqr),
        parameter_correlation_matrix: generic.parameter_correlation_matrix,
        iterations: generic.iterations,
        converged: generic.converged,
        weight: generic.weight,
    })
}

pub fn complete_initial_params(
    data: &EisData,
    partial: PartialRqrInit,
    auto: bool,
) -> Result<RqrParams> {
    let model: EcmModelSpec = "R_QR".parse()?;
    let params = complete_initial_ecm(
        data,
        &model,
        PartialEcmInit {
            rs: partial.rs,
            r1: partial.rct,
            q1: partial.q,
            n1: partial.n,
            ..PartialEcmInit::default()
        },
        auto,
    )?;
    Ok(ecm_to_rqr(&params))
}

pub fn auto_init_rqr(data: &EisData) -> RqrParams {
    let model: EcmModelSpec = "R_QR".parse().expect("built-in model must parse");
    ecm_to_rqr(&auto_init_ecm(data, &model))
}

fn rqr_to_ecm(params: RqrParams) -> EcmParams {
    EcmParams {
        rs: params.rs,
        branches: vec![ParallelBranchParams {
            r: params.rct,
            reactive: ReactiveElementParams::Cpe {
                q: params.q,
                n: params.n,
            },
        }],
        warburg_sigma: None,
    }
}

fn ecm_to_rqr(params: &EcmParams) -> RqrParams {
    let ReactiveElementParams::Cpe { q, n } = params.branches[0].reactive else {
        unreachable!("R_QR must contain a CPE branch")
    };
    RqrParams {
        rs: params.rs,
        rct: params.branches[0].r,
        q,
        n,
    }
}

fn vector_to_rqr(values: &[f64]) -> RqrParams {
    RqrParams {
        rs: values[0],
        rct: values[1],
        q: values[2],
        n: values[3],
    }
}
