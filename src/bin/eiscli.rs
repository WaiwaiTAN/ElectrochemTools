use ElectrochemTools::drt::{DrtSettings, TauGridMode, scan_lambda, solve_drt};
use ElectrochemTools::drt_compare::compare_with_matlab_outputs;
use ElectrochemTools::ecm::RqrParams;
use ElectrochemTools::eis_io::{read_eis_with_cleaning, write_impedance_csv};
use ElectrochemTools::fit::{
    PartialRqrInit, RqrFitSettings, Weighting, complete_initial_params, fit_rqr,
};
use ElectrochemTools::plot::{write_drt_gamma_svg, write_nyquist_svg};
use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use serde::Serialize;
use std::fs;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "EIS post-processing CLI for DRT and ECM fitting"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Tikhonov DRT MVP using direct Debye discretization.
    Drt {
        input: PathBuf,
        #[arg(long, default_value_t = 1.0e-3)]
        lambda: f64,
        #[arg(long)]
        auto_lambda: bool,
        #[arg(long, default_value_t = 1.0e-6)]
        lambda_min: f64,
        #[arg(long, default_value_t = 1.0)]
        lambda_max: f64,
        #[arg(long, default_value_t = 20)]
        n_lambda: usize,
        #[arg(long)]
        tau_min: Option<f64>,
        #[arg(long)]
        tau_max: Option<f64>,
        #[arg(long, default_value_t = 100)]
        n_tau: usize,
        #[arg(long, value_enum, default_value_t = TauGridMode::Logspace)]
        tau_grid: TauGridMode,
        #[arg(long, default_value_t = 1)]
        regularization_order: usize,
        #[arg(long)]
        flip_imag: bool,
        #[arg(long)]
        drop_positive_imag: bool,
        #[arg(long)]
        nonnegative: bool,
        #[arg(long)]
        fit_inductance: bool,
        #[arg(long)]
        credible_intervals: bool,
        #[arg(long)]
        compare_matlab_drt: Option<PathBuf>,
        #[arg(long)]
        compare_matlab_regression: Option<PathBuf>,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Equivalent circuit fitting. Currently supports R_QR.
    FitEcm {
        input: PathBuf,
        #[arg(long)]
        model: String,
        #[arg(long)]
        rs: Option<f64>,
        #[arg(long)]
        rct: Option<f64>,
        #[arg(long)]
        q: Option<f64>,
        #[arg(long)]
        n: Option<f64>,
        #[arg(long)]
        auto_init: bool,
        #[arg(long, value_enum, default_value_t = Weighting::Proportional)]
        weight: Weighting,
        #[arg(long, default_value_t = 200)]
        max_iter: usize,
        #[arg(long, default_value_t = 1.0e-10)]
        tol: f64,
        #[arg(long)]
        flip_imag: bool,
        #[arg(long)]
        drop_positive_imag: bool,
        #[arg(long)]
        include_correlation_matrix: bool,
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

#[derive(Serialize)]
struct DrtSummary {
    lambda: f64,
    tau_min: f64,
    tau_max: f64,
    n_tau: usize,
    tau_grid: TauGridMode,
    fit_inductance: bool,
    regularization_order: usize,
    nonnegative: bool,
    n_points: usize,
    r_inf: f64,
    inductance: f64,
    polarization_resistance: f64,
    rmse_real: f64,
    rmse_imag: f64,
    rmse_magnitude: f64,
    relative_rmse_percent: f64,
    kk_real_to_imag_relative_rmse_percent: f64,
    kk_imag_to_real_relative_rmse_percent: f64,
    kk_mean_score: f64,
    credible_intervals: bool,
    inductance_std: Option<f64>,
    r_inf_std: Option<f64>,
    note: String,
}

#[derive(Serialize)]
struct FitParamsSummary {
    model: String,
    rs: f64,
    rct: f64,
    q: f64,
    n: f64,
    weighted_sse: f64,
    mean_weighted_chi_square: f64,
    reduced_chi_square: f64,
    relative_rmse_percent: f64,
    rmse_real: f64,
    rmse_imag: f64,
    rmse_magnitude: f64,
    parameter_std_errors: Option<RqrParams>,
    parameter_rel_std_error_percent: Option<RqrParams>,
    parameter_correlation_labels: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    parameter_correlation_matrix: Option<Vec<Vec<f64>>>,
    weight: Weighting,
    n_points: usize,
    iterations: usize,
    converged: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Drt {
            input,
            lambda,
            auto_lambda,
            lambda_min,
            lambda_max,
            n_lambda,
            tau_min,
            tau_max,
            n_tau,
            tau_grid,
            regularization_order,
            flip_imag,
            drop_positive_imag,
            nonnegative,
            fit_inductance,
            credible_intervals,
            compare_matlab_drt,
            compare_matlab_regression,
            out,
        } => run_drt(
            input,
            lambda,
            auto_lambda,
            lambda_min,
            lambda_max,
            n_lambda,
            tau_min,
            tau_max,
            n_tau,
            tau_grid,
            regularization_order,
            flip_imag,
            drop_positive_imag,
            nonnegative,
            fit_inductance,
            credible_intervals,
            compare_matlab_drt,
            compare_matlab_regression,
            out,
        ),
        Commands::FitEcm {
            input,
            model,
            rs,
            rct,
            q,
            n,
            auto_init,
            weight,
            max_iter,
            tol,
            flip_imag,
            drop_positive_imag,
            include_correlation_matrix,
            out,
        } => run_fit_ecm(
            input,
            model,
            rs,
            rct,
            q,
            n,
            auto_init,
            weight,
            max_iter,
            tol,
            flip_imag,
            drop_positive_imag,
            include_correlation_matrix,
            out,
        ),
    }
}

fn run_drt(
    input: PathBuf,
    lambda: f64,
    auto_lambda: bool,
    lambda_min: f64,
    lambda_max: f64,
    n_lambda: usize,
    tau_min: Option<f64>,
    tau_max: Option<f64>,
    n_tau: usize,
    tau_grid: TauGridMode,
    regularization_order: usize,
    flip_imag: bool,
    drop_positive_imag: bool,
    nonnegative: bool,
    fit_inductance: bool,
    credible_intervals: bool,
    compare_matlab_drt: Option<PathBuf>,
    compare_matlab_regression: Option<PathBuf>,
    out: Option<PathBuf>,
) -> Result<()> {
    let mut data = read_eis_with_cleaning(&input, drop_positive_imag)?;
    if flip_imag {
        data.flip_imag();
    }
    data.warn_if_imag_mostly_positive();
    let out = out.unwrap_or_else(|| default_output_dir(&input, "drt"));
    fs::create_dir_all(&out).with_context(|| format!("failed to create {}", out.display()))?;

    let mut settings = DrtSettings {
        lambda,
        tau_min,
        tau_max,
        n_tau,
        tau_grid,
        fit_inductance,
        regularization_order,
        nonnegative,
        credible_intervals,
    };

    if auto_lambda {
        let (best_lambda, scan) = scan_lambda(&data, &settings, lambda_min, lambda_max, n_lambda)?;
        settings.lambda = best_lambda;
        let mut scan_writer = csv::Writer::from_path(out.join("lambda_scan.csv"))?;
        scan_writer.write_record(["lambda", "relative_rmse", "roughness", "score"])?;
        for point in scan {
            scan_writer.serialize(point)?;
        }
        scan_writer.flush()?;
    }

    let result = solve_drt(&data, &settings)?;

    let mut gamma_writer = csv::Writer::from_path(out.join("gamma.csv"))?;
    gamma_writer.write_record(["tau", "log10_tau", "gamma"])?;
    for (&tau, &gamma) in result.tau.iter().zip(&result.gamma) {
        gamma_writer.serialize((tau, tau.log10(), gamma))?;
    }
    gamma_writer.flush()?;

    if let Some(ci) = &result.credible_intervals {
        let mut ci_writer = csv::Writer::from_path(out.join("gamma_ci.csv"))?;
        ci_writer.write_record([
            "tau",
            "log10_tau",
            "gamma",
            "gamma_std",
            "gamma_lower_95",
            "gamma_upper_95",
        ])?;
        for idx in 0..result.tau.len() {
            ci_writer.serialize((
                result.tau[idx],
                result.tau[idx].log10(),
                result.gamma[idx],
                ci.gamma_std[idx],
                ci.gamma_lower_95[idx],
                ci.gamma_upper_95[idx],
            ))?;
        }
        ci_writer.flush()?;
    }

    let mut peak_writer = csv::Writer::from_path(out.join("drt_peaks.csv"))?;
    peak_writer.write_record(["tau", "log10_tau", "gamma"])?;
    for peak in &result.peaks {
        peak_writer.serialize(peak)?;
    }
    peak_writer.flush()?;

    let mut drttools_writer = csv::WriterBuilder::new()
        .has_headers(false)
        .from_path(out.join("drttools_compatible_drt.csv"))?;
    drttools_writer.write_record(["L", &format_matlab_exp(result.inductance)])?;
    drttools_writer.write_record(["R", &format_matlab_exp(result.polarization_resistance)])?;
    drttools_writer.write_record(["tau", "gamma(tau)"])?;
    for (&tau, &gamma) in result.tau.iter().zip(&result.gamma) {
        drttools_writer.write_record([format_matlab_exp(tau), format_matlab_exp(gamma)])?;
    }
    drttools_writer.flush()?;

    write_impedance_csv(
        &out.join("reconstructed_impedance.csv"),
        &data,
        &result.z_fit,
    )?;
    write_drt_svgs(&out, &data, &result)?;
    let mut kk_writer = csv::Writer::from_path(out.join("kk_consistency.csv"))?;
    kk_writer.write_record([
        "frequency",
        "Z_real_exp",
        "Z_imag_exp",
        "Z_real_from_imag",
        "Z_imag_from_real",
        "residual_real_from_imag",
        "residual_imag_from_real",
    ])?;
    for row in &result.kk.rows {
        kk_writer.serialize(row)?;
    }
    kk_writer.flush()?;
    fs::write(
        out.join("kk_summary.json"),
        serde_json::to_string_pretty(&result.kk)?,
    )?;
    let summary = DrtSummary {
        lambda: result.settings_used.lambda,
        tau_min: result.settings_used.tau_min,
        tau_max: result.settings_used.tau_max,
        n_tau: result.settings_used.n_tau,
        tau_grid: result.settings_used.tau_grid,
        fit_inductance: result.settings_used.fit_inductance,
        regularization_order: result.settings_used.regularization_order,
        nonnegative: result.settings_used.nonnegative,
        n_points: data.len(),
        r_inf: result.r_inf,
        inductance: result.inductance,
        polarization_resistance: result.polarization_resistance,
        rmse_real: result.metrics.rmse_real,
        rmse_imag: result.metrics.rmse_imag,
        rmse_magnitude: result.metrics.rmse_magnitude,
        relative_rmse_percent: result.metrics.relative_rmse_percent,
        kk_real_to_imag_relative_rmse_percent: result
            .kk
            .real_to_imag_relative_rmse_percent,
        kk_imag_to_real_relative_rmse_percent: result
            .kk
            .imag_to_real_relative_rmse_percent,
        kk_mean_score: result.kk.mean_score,
        credible_intervals: result.settings_used.credible_intervals,
        inductance_std: result
            .credible_intervals
            .as_ref()
            .and_then(|intervals| intervals.inductance_std),
        r_inf_std: result
            .credible_intervals
            .as_ref()
            .map(|intervals| intervals.r_inf_std),
        note: if result.settings_used.nonnegative {
            "bounded active-set nonnegative Tikhonov DRT; credible intervals are a linear-Gaussian approximation when requested"
        } else {
            "unconstrained Tikhonov DRT; use --nonnegative to enforce bounded nonnegative coefficients"
        }
        .to_string(),
    };
    fs::write(
        out.join("residual_summary.json"),
        serde_json::to_string_pretty(&summary)?,
    )?;
    if compare_matlab_drt.is_some() || compare_matlab_regression.is_some() {
        let comparison = compare_with_matlab_outputs(
            &data,
            &result,
            compare_matlab_drt.as_deref(),
            compare_matlab_regression.as_deref(),
        )?;
        fs::write(
            out.join("matlab_comparison.json"),
            serde_json::to_string_pretty(&comparison)?,
        )?;
    }
    println!("DRT analysis complete: {}", out.display());
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_fit_ecm(
    input: PathBuf,
    model: String,
    rs: Option<f64>,
    rct: Option<f64>,
    q: Option<f64>,
    n: Option<f64>,
    auto_init: bool,
    weight: Weighting,
    max_iter: usize,
    tol: f64,
    flip_imag: bool,
    drop_positive_imag: bool,
    include_correlation_matrix: bool,
    out: Option<PathBuf>,
) -> Result<()> {
    if model.to_ascii_uppercase() != "R_QR" {
        bail!("unsupported model '{}'; currently supported: R_QR", model);
    }

    let mut data = read_eis_with_cleaning(&input, drop_positive_imag)?;
    if flip_imag {
        data.flip_imag();
    }
    data.warn_if_imag_mostly_positive();
    let out = out.unwrap_or_else(|| default_output_dir(&input, "ecm"));
    fs::create_dir_all(&out).with_context(|| format!("failed to create {}", out.display()))?;

    let initial: RqrParams =
        complete_initial_params(&data, PartialRqrInit { rs, rct, q, n }, auto_init)?;
    let result = fit_rqr(
        &data,
        &RqrFitSettings {
            initial,
            weight,
            max_iter,
            tol,
        },
    )?;

    write_impedance_csv(&out.join("fitted_impedance.csv"), &data, &result.z_fit)?;
    write_nyquist_svg(
        &out.join("nyquist_fit.svg"),
        "R(QR) Nyquist Fit",
        &data.z_real,
        &data.z_imag,
        &result.z_fit,
    )?;
    let summary = FitParamsSummary {
        model: "R_QR".to_string(),
        rs: result.params.rs,
        rct: result.params.rct,
        q: result.params.q,
        n: result.params.n,
        weighted_sse: result.weighted_sse,
        mean_weighted_chi_square: result.mean_weighted_chi_square,
        reduced_chi_square: result.reduced_chi_square,
        relative_rmse_percent: result.metrics.relative_rmse_percent,
        rmse_real: result.metrics.rmse_real,
        rmse_imag: result.metrics.rmse_imag,
        rmse_magnitude: result.metrics.rmse_magnitude,
        parameter_std_errors: result.parameter_std_errors,
        parameter_rel_std_error_percent: result.parameter_rel_std_error_percent,
        parameter_correlation_labels: vec![
            "Rs".to_string(),
            "R".to_string(),
            "Q".to_string(),
            "alpha".to_string(),
        ],
        parameter_correlation_matrix: include_correlation_matrix
            .then_some(result.parameter_correlation_matrix.clone())
            .flatten(),
        weight: result.weight,
        n_points: data.len(),
        iterations: result.iterations,
        converged: result.converged,
    };
    fs::write(
        out.join("fit_params.json"),
        serde_json::to_string_pretty(&summary)?,
    )?;
    print_fit_summary(&result);
    println!("ECM fitting complete: {}", out.display());
    Ok(())
}

fn write_drt_svgs(
    out: &PathBuf,
    data: &ElectrochemTools::types::EisData,
    result: &ElectrochemTools::drt::DrtResult,
) -> Result<()> {
    write_drt_gamma_svg(
        &out.join("drt_gamma.svg"),
        &result.tau,
        &result.gamma,
        "DRT Gamma",
    )?;
    write_nyquist_svg(
        &out.join("nyquist_reconstruction.svg"),
        "DRT Nyquist Reconstruction",
        &data.z_real,
        &data.z_imag,
        &result.z_fit,
    )
}

fn print_fit_summary(result: &ElectrochemTools::fit::RqrFitResult) {
    let rel = result.parameter_rel_std_error_percent;
    println!("Fit quality:");
    println!("  weighted SSE              = {:.3e}", result.weighted_sse);
    println!(
        "  mean weighted chi-square  = {:.3e}",
        result.mean_weighted_chi_square
    );
    println!(
        "  reduced chi-square        = {:.3e}",
        result.reduced_chi_square
    );
    println!(
        "  relative RMSE             = {:.2} %",
        result.metrics.relative_rmse_percent
    );
    println!();
    println!("Parameters:");
    print_param_line("Rs", result.params.rs, "ohm", rel.map(|value| value.rs));
    print_param_line("R", result.params.rct, "ohm", rel.map(|value| value.rct));
    print_param_line("Q", result.params.q, "", rel.map(|value| value.q));
    print_param_line("alpha", result.params.n, "", rel.map(|value| value.n));
}

fn default_output_dir(input: &std::path::Path, suffix: &str) -> PathBuf {
    let stem = input
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("eis");
    let dir_name = format!("{stem}_{suffix}");
    input
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| std::path::Path::new("."))
        .join(dir_name)
}

fn print_param_line(name: &str, value: f64, unit: &str, rel_error: Option<f64>) {
    let unit_text = if unit.is_empty() {
        String::new()
    } else {
        format!(" {unit}")
    };
    let rel_text = rel_error
        .map(|value| format!("{value:.1} %"))
        .unwrap_or_else(|| "n/a".to_string());
    println!(
        "  {name:<7} = {:<12}{unit_text:<8} rel. std. error = {rel_text}",
        format_compact(value),
    );
}

fn format_compact(value: f64) -> String {
    if value == 0.0 || (value.abs() >= 1.0e-3 && value.abs() < 1.0e4) {
        format!("{value:.5}")
    } else {
        format!("{value:.3e}")
    }
}

fn format_matlab_exp(value: f64) -> String {
    let raw = format!("{value:.6e}");
    let Some((mantissa, exponent)) = raw.split_once('e') else {
        return raw;
    };
    let exp_value = exponent.parse::<i32>().unwrap_or(0);
    format!("{mantissa}e{exp_value:+03}")
}
