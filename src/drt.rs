mod bayesian;
pub mod discretization;
pub mod kernel;
pub mod regularization;
pub mod solver;

pub use crate::regularization::SolverReport;
use crate::types::{EisData, FitMetrics, calculate_fit_metrics};
use anyhow::{Result, bail};
pub use bayesian::BayesianSettings;
use bayesian::sample_nonnegative_posterior;
pub use discretization::{DrtBasis, ShapeControl};
use discretization::{
    DrtDiscretization, GaussianDiscretization, PiecewiseLinearDiscretization,
    evaluate_gaussian_profile,
};
use kernel::{
    BasisKernelMatrices, assemble_combined_from_kernels, assemble_imag_from_kernels,
    assemble_real_from_kernels, reconstruct_from_kernels,
};
pub use kernel::{delta_ln_tau, reconstruct_impedance, reconstruct_impedance_with_inductance};
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;
use serde::Serialize;
use solver::solve_coefficients;
pub use solver::{DrtConstraintConfig, DrtSolverOptions};
use std::f64::consts::PI;

const DRTTOOLS_PLOT_MARGIN_DECADES: f64 = 0.5;

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
    pub basis: DrtBasis,
    pub shape_control: ShapeControl,
    pub shape_coefficient: f64,
    pub fit_inductance: bool,
    pub regularization_order: usize,
    pub nonnegative: bool,
    pub credible_intervals: bool,
    pub bayesian: Option<BayesianSettings>,
    pub solver: DrtSolverOptions,
}

#[derive(Debug, Clone)]
pub struct DrtResult {
    pub tau: Vec<f64>,
    pub gamma: Vec<f64>,
    /// Tau coordinates used only for smooth visualization of the fitted profile.
    pub plot_tau: Vec<f64>,
    /// Gamma evaluated at `plot_tau`; solver-grid outputs remain in `tau` and `gamma`.
    pub plot_gamma: Vec<f64>,
    pub r_inf: f64,
    pub z_fit: Vec<Complex<f64>>,
    pub settings_used: DrtSettingsUsed,
    pub metrics: FitMetrics,
    pub polarization_resistance: f64,
    pub inductance: f64,
    pub peaks: Vec<DrtPeak>,
    pub credible_intervals: Option<DrtCredibleIntervals>,
    pub bayesian: Option<DrtBayesianResult>,
    pub kk: KkConsistencyResult,
    pub solver_report: SolverReport,
}

#[derive(Debug, Clone, Serialize)]
pub struct DrtSettingsUsed {
    pub lambda: f64,
    pub tau_min: f64,
    pub tau_max: f64,
    pub n_tau: usize,
    pub tau_grid: TauGridMode,
    pub basis: DrtBasis,
    pub shape_control: ShapeControl,
    pub shape_coefficient: f64,
    pub epsilon: Option<f64>,
    pub fit_inductance: bool,
    pub regularization_order: usize,
    pub nonnegative: bool,
    pub credible_intervals: bool,
    pub solver_max_iterations: usize,
    pub solver_tolerance: f64,
    pub constraints: DrtConstraintConfig,
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
pub struct DrtBayesianResult {
    pub chains: usize,
    pub samples_per_chain: usize,
    pub burn_in: usize,
    pub retained_samples_per_chain: usize,
    pub total_samples: usize,
    pub total_retained_samples: usize,
    pub seed: u64,
    pub chain_seeds: Vec<u64>,
    pub lower_probability: f64,
    pub upper_probability: f64,
    pub noise_std: f64,
    pub bounce_count: usize,
    pub coefficient_r_hat: Vec<f64>,
    pub coefficient_effective_sample_size: Vec<f64>,
    pub max_r_hat: f64,
    pub min_effective_sample_size: f64,
    pub diagnostics_qualified: bool,
    pub coefficient_mean: Vec<f64>,
    pub coefficient_lower: Vec<f64>,
    pub coefficient_upper: Vec<f64>,
    pub gamma_mean: Vec<f64>,
    pub gamma_lower: Vec<f64>,
    pub gamma_upper: Vec<f64>,
    #[serde(skip_serializing)]
    pub plot_gamma_mean: Vec<f64>,
    #[serde(skip_serializing)]
    pub plot_gamma_lower: Vec<f64>,
    #[serde(skip_serializing)]
    pub plot_gamma_upper: Vec<f64>,
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
    let (tau, used_tau_grid) = match settings.basis {
        DrtBasis::PiecewiseLinear => (
            match settings.tau_grid {
                TauGridMode::Logspace => make_log_tau_grid(tau_min, tau_max, settings.n_tau)?,
                TauGridMode::Drttools => make_drttools_tau_grid(data)?,
            },
            settings.tau_grid,
        ),
        DrtBasis::Gaussian => (make_drttools_tau_grid(data)?, TauGridMode::Drttools),
    };
    let discretization: Box<dyn DrtDiscretization> = match settings.basis {
        DrtBasis::PiecewiseLinear => Box::new(PiecewiseLinearDiscretization::new(tau)?),
        DrtBasis::Gaussian => Box::new(GaussianDiscretization::new(
            tau,
            settings.shape_control,
            settings.shape_coefficient,
        )?),
    };
    let tau = discretization.tau();
    let n_tau = tau.len();
    let n_unpenalized = if settings.fit_inductance { 2 } else { 1 };
    let gamma_offset = n_unpenalized;
    let kernels = discretization.kernel_matrices(&data.frequency_hz)?;
    let (a, b) = assemble_combined_from_kernels(data, &kernels, settings.fit_inductance);
    let penalty = discretization.penalty(settings.regularization_order, n_unpenalized)?;
    let solution = solve_coefficients(
        &a,
        &b,
        settings.lambda,
        &penalty,
        settings.nonnegative,
        settings.fit_inductance,
        settings.solver,
    )?;
    let x = solution.coefficients;
    let (inductance, r_inf) = if settings.fit_inductance {
        (x[0], x[1])
    } else {
        (0.0, x[0])
    };
    let basis_coefficients = DVector::from_iterator(n_tau, x.iter().skip(gamma_offset).copied());
    let gamma_vector = discretization.map_coefficients_to_gamma(&basis_coefficients)?;
    let gamma: Vec<f64> = gamma_vector.iter().copied().collect();
    let basis_coefficients_vec = basis_coefficients.iter().copied().collect::<Vec<_>>();
    let (plot_tau, plot_gamma) = match discretization.basis() {
        DrtBasis::PiecewiseLinear => (tau.to_vec(), gamma.clone()),
        DrtBasis::Gaussian => {
            let plot_tau_min = 10_f64.powf(tau[0].log10() - DRTTOOLS_PLOT_MARGIN_DECADES);
            let plot_tau_max = 10_f64.powf(tau[n_tau - 1].log10() + DRTTOOLS_PLOT_MARGIN_DECADES);
            let plot_tau = make_log_tau_grid(plot_tau_min, plot_tau_max, 10 * n_tau)?;
            let plot_gamma = evaluate_gaussian_profile(
                tau,
                &basis_coefficients_vec,
                discretization
                    .epsilon()
                    .expect("Gaussian discretization must provide epsilon"),
                &plot_tau,
            )?;
            (plot_tau, plot_gamma)
        }
    };
    let z_fit = match discretization.basis() {
        DrtBasis::PiecewiseLinear => reconstruct_impedance_with_inductance(
            &data.frequency_hz,
            r_inf,
            inductance,
            tau,
            &gamma,
        ),
        DrtBasis::Gaussian => reconstruct_from_kernels(
            &data.frequency_hz,
            r_inf,
            inductance,
            &kernels,
            &basis_coefficients_vec,
        ),
    };
    let metrics = calculate_fit_metrics(data, &z_fit)?;
    let polarization_resistance = polarization_resistance(tau, &gamma);
    let peaks = detect_drt_peaks(tau, &gamma, 0.01);
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
            &discretization.gamma_mapping_matrix(),
            &gamma,
        )
    } else {
        None
    };
    let bayesian = if let Some(bayesian_settings) = settings.bayesian {
        let sampling = sample_nonnegative_posterior(
            &a,
            &b,
            &x,
            settings.lambda,
            &penalty,
            n_unpenalized,
            bayesian_settings,
        )?;
        let mapping = discretization.gamma_mapping_matrix();
        let gamma_mean = &mapping * DVector::from_vec(sampling.coefficient_mean.clone());
        let gamma_lower = &mapping * DVector::from_vec(sampling.coefficient_lower.clone());
        let gamma_upper = &mapping * DVector::from_vec(sampling.coefficient_upper.clone());
        let (plot_gamma_mean, plot_gamma_lower, plot_gamma_upper) = match discretization.basis() {
            DrtBasis::PiecewiseLinear => (
                gamma_mean.iter().copied().collect(),
                gamma_lower.iter().copied().collect(),
                gamma_upper.iter().copied().collect(),
            ),
            DrtBasis::Gaussian => {
                let epsilon = discretization
                    .epsilon()
                    .expect("Gaussian discretization must provide epsilon");
                (
                    evaluate_gaussian_profile(tau, &sampling.coefficient_mean, epsilon, &plot_tau)?,
                    evaluate_gaussian_profile(
                        tau,
                        &sampling.coefficient_lower,
                        epsilon,
                        &plot_tau,
                    )?,
                    evaluate_gaussian_profile(
                        tau,
                        &sampling.coefficient_upper,
                        epsilon,
                        &plot_tau,
                    )?,
                )
            }
        };
        Some(DrtBayesianResult {
            chains: sampling.chains,
            samples_per_chain: sampling.samples_per_chain,
            burn_in: sampling.burn_in,
            retained_samples_per_chain: sampling.retained_samples_per_chain,
            total_samples: sampling.total_samples,
            total_retained_samples: sampling.total_retained_samples,
            seed: sampling.seed,
            chain_seeds: sampling.chain_seeds,
            lower_probability: sampling.lower_probability,
            upper_probability: sampling.upper_probability,
            noise_std: sampling.noise_std,
            bounce_count: sampling.bounce_count,
            coefficient_r_hat: sampling.coefficient_r_hat,
            coefficient_effective_sample_size: sampling.coefficient_effective_sample_size,
            max_r_hat: sampling.max_r_hat,
            min_effective_sample_size: sampling.min_effective_sample_size,
            diagnostics_qualified: sampling.diagnostics_qualified,
            coefficient_mean: sampling.coefficient_mean,
            coefficient_lower: sampling.coefficient_lower,
            coefficient_upper: sampling.coefficient_upper,
            gamma_mean: gamma_mean.iter().copied().collect(),
            gamma_lower: gamma_lower.iter().copied().collect(),
            gamma_upper: gamma_upper.iter().copied().collect(),
            plot_gamma_mean,
            plot_gamma_lower,
            plot_gamma_upper,
            note: "DRTtools-style exact HMC samples of the Gaussian posterior truncated to nonnegative DRT basis coefficients; bounds are R-5 0.5% and 99.5% coefficient quantiles mapped to gamma"
                .to_string(),
        })
    } else {
        None
    };
    let kk = estimate_kk_consistency_with_discretization(
        data,
        discretization.as_ref(),
        &kernels,
        settings.lambda,
        settings.regularization_order,
    )?;

    Ok(DrtResult {
        tau: tau.to_vec(),
        gamma,
        plot_tau,
        plot_gamma,
        r_inf,
        z_fit,
        settings_used: DrtSettingsUsed {
            lambda: settings.lambda,
            tau_min: used_tau_min,
            tau_max: used_tau_max,
            n_tau,
            tau_grid: used_tau_grid,
            basis: discretization.basis(),
            shape_control: settings.shape_control,
            shape_coefficient: settings.shape_coefficient,
            epsilon: discretization.epsilon(),
            fit_inductance: settings.fit_inductance,
            regularization_order: settings.regularization_order,
            nonnegative: settings.nonnegative,
            credible_intervals: settings.credible_intervals,
            solver_max_iterations: settings.solver.max_iterations,
            solver_tolerance: settings.solver.tolerance,
            constraints: settings.solver.constraints,
        },
        metrics,
        polarization_resistance,
        inductance,
        peaks,
        credible_intervals,
        bayesian,
        kk,
        solver_report: solution.report,
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
        trial.bayesian = None;
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
    gamma_mapping: &DMatrix<f64>,
    gamma: &[f64],
) -> Option<DrtCredibleIntervals> {
    let residual = a * x - b;
    let dof = (b.len() as f64 - x.len() as f64).max(1.0);
    let sigma2 = residual.dot(&residual) / dof;
    let precision = a.transpose() * a + penalty.scale(lambda);
    let covariance = precision.try_inverse()?.scale(sigma2);
    let coefficient_covariance = covariance
        .view((n_unpenalized, n_unpenalized), (n_gamma, n_gamma))
        .into_owned();
    let gamma_covariance = gamma_mapping * coefficient_covariance * gamma_mapping.transpose();
    let mut gamma_std = Vec::with_capacity(n_gamma);
    let mut gamma_lower_95 = Vec::with_capacity(n_gamma);
    let mut gamma_upper_95 = Vec::with_capacity(n_gamma);
    for idx in 0..n_gamma {
        let std = gamma_covariance[(idx, idx)].max(0.0).sqrt();
        let center = gamma[idx];
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
    let discretization = PiecewiseLinearDiscretization::new(tau.to_vec())?;
    let kernels = discretization.kernel_matrices(&data.frequency_hz)?;
    estimate_kk_consistency_with_discretization(data, &discretization, &kernels, lambda, order)
}

fn estimate_kk_consistency_with_discretization(
    data: &EisData,
    discretization: &dyn DrtDiscretization,
    kernels: &BasisKernelMatrices,
    lambda: f64,
    order: usize,
) -> Result<KkConsistencyResult> {
    let (a_re, b_re) = assemble_real_from_kernels(data, kernels);
    let (a_im, b_im) = assemble_imag_from_kernels(data, kernels);
    let penalty = discretization.penalty(order, 1)?;
    let x_re = solve_coefficients(
        &a_re,
        &b_re,
        lambda,
        &penalty,
        false,
        false,
        DrtSolverOptions::default(),
    )?
    .coefficients;
    let x_im = solve_coefficients(
        &a_im,
        &b_im,
        lambda,
        &penalty,
        false,
        false,
        DrtSolverOptions::default(),
    )?
    .coefficients;
    let coefficients_re: Vec<f64> = x_re.iter().skip(1).copied().collect();
    let coefficients_im: Vec<f64> = x_im.iter().skip(1).copied().collect();

    let (z_from_re, z_real_no_offset) = match discretization.basis() {
        DrtBasis::PiecewiseLinear => (
            reconstruct_impedance(
                &data.frequency_hz,
                x_re[0],
                discretization.tau(),
                &coefficients_re,
            ),
            reconstruct_impedance(
                &data.frequency_hz,
                0.0,
                discretization.tau(),
                &coefficients_im,
            ),
        ),
        DrtBasis::Gaussian => (
            reconstruct_from_kernels(&data.frequency_hz, x_re[0], 0.0, kernels, &coefficients_re),
            reconstruct_from_kernels(&data.frequency_hz, 0.0, 0.0, kernels, &coefficients_im),
        ),
    };
    let r_inf_from_imag = data
        .z_real
        .iter()
        .zip(&z_real_no_offset)
        .map(|(&exp, fit)| exp - fit.re)
        .sum::<f64>()
        / data.len() as f64;
    let z_from_im = match discretization.basis() {
        DrtBasis::PiecewiseLinear => reconstruct_impedance(
            &data.frequency_hz,
            r_inf_from_imag,
            discretization.tau(),
            &coefficients_im,
        ),
        DrtBasis::Gaussian => reconstruct_from_kernels(
            &data.frequency_hz,
            r_inf_from_imag,
            0.0,
            kernels,
            &coefficients_im,
        ),
    };

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
