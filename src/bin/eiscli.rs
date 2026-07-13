use anyhow::{Context, Result, bail};
use clap::{Args, Parser, Subcommand};
use electrochem_tools::batch::{
    BatchOptions, BatchReport, BatchStatus, default_jobs, run_batch, run_batch_with_resume,
};
use electrochem_tools::drt::{
    DrtBasis, DrtConstraintConfig, DrtSettings, DrtSolverOptions, ShapeControl, SolverReport,
    TauGridMode, scan_lambda, solve_drt,
};
use electrochem_tools::drt_compare::compare_with_matlab_outputs;
use electrochem_tools::ecm::{EcmModelSpec, EcmParams};
use electrochem_tools::eis::{CleanOptions, ImagSignPolicy, clean_file_to};
use electrochem_tools::eis_io::{read_eis_with_cleaning, write_impedance_csv};
use electrochem_tools::fit::{
    EcmFitSettings, PartialEcmInit, Weighting, complete_initial_ecm, fit_ecm,
};
use electrochem_tools::plot::{write_drt_gamma_svg, write_nyquist_svg};
use electrochem_tools::run_manifest::{
    RunManifest, collect_output_files, execute_manifested, verify_resume,
};
use electrochem_tools::types::EisData;
use serde::Serialize;
use serde_json::json;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

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

#[derive(Args, Debug, Clone)]
struct BatchArgs {
    #[arg(long)]
    jobs: Option<usize>,
    #[arg(long)]
    fail_fast: bool,
    #[arg(long)]
    overwrite: bool,
    #[arg(long)]
    resume: bool,
    #[arg(long)]
    out_root: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
struct CleanBatchArgs {
    #[arg(long)]
    jobs: Option<usize>,
    #[arg(long)]
    fail_fast: bool,
    #[arg(long)]
    overwrite: bool,
    #[arg(long)]
    out_root: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
struct EcmInitialArgs {
    #[arg(long)]
    rs: Option<f64>,
    #[arg(long, visible_alias = "rct")]
    r1: Option<f64>,
    #[arg(long, visible_alias = "c")]
    c1: Option<f64>,
    #[arg(long, visible_alias = "q")]
    q1: Option<f64>,
    #[arg(long, visible_alias = "n")]
    n1: Option<f64>,
    #[arg(long)]
    r2: Option<f64>,
    #[arg(long)]
    c2: Option<f64>,
    #[arg(long)]
    q2: Option<f64>,
    #[arg(long)]
    n2: Option<f64>,
    #[arg(long, visible_alias = "sigma-w")]
    warburg_sigma: Option<f64>,
}

impl From<EcmInitialArgs> for PartialEcmInit {
    fn from(value: EcmInitialArgs) -> Self {
        Self {
            rs: value.rs,
            r1: value.r1,
            c1: value.c1,
            q1: value.q1,
            n1: value.n1,
            r2: value.r2,
            c2: value.c2,
            q2: value.q2,
            n2: value.n2,
            warburg_sigma: value.warburg_sigma,
        }
    }
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Validate and clean one or more EIS files using the shared input layer.
    Clean {
        #[arg(short = 'i', long = "input", required = true, num_args = 1..)]
        input: Vec<PathBuf>,
        #[arg(long)]
        lenient: bool,
        #[arg(long, value_enum, default_value_t = ImagSignPolicy::Preserve)]
        imag_sign: ImagSignPolicy,
        #[arg(long)]
        keep_positive_imag: bool,
        #[command(flatten)]
        batch: CleanBatchArgs,
    },
    /// Tikhonov DRT with piecewise-linear or Gaussian discretization.
    Drt {
        #[arg(short = 'i', long = "input", required = true, num_args = 1..)]
        input: Vec<PathBuf>,
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
        #[arg(long, value_enum, default_value_t = DrtBasis::PiecewiseLinear)]
        basis: DrtBasis,
        #[arg(long, value_enum, default_value_t = ShapeControl::Fwhm)]
        shape_control: ShapeControl,
        #[arg(long, default_value_t = 0.5)]
        shape_coefficient: f64,
        #[arg(long, default_value_t = 1)]
        regularization_order: usize,
        #[arg(long)]
        flip_imag: bool,
        #[arg(long)]
        keep_positive_imag: bool,
        #[arg(long)]
        nonnegative: bool,
        #[arg(long)]
        fit_inductance: bool,
        #[arg(long)]
        credible_intervals: bool,
        #[arg(long, default_value_t = 1_000)]
        solver_max_iterations: usize,
        #[arg(long, default_value_t = 1.0e-9)]
        solver_tolerance: f64,
        #[arg(long)]
        allow_negative_r_inf: bool,
        #[arg(long)]
        nonnegative_inductance: bool,
        #[arg(long)]
        compare_matlab_drt: Option<PathBuf>,
        #[arg(long)]
        compare_matlab_regression: Option<PathBuf>,
        #[command(flatten)]
        batch: BatchArgs,
    },
    /// Equivalent-circuit fitting with one/two RC or RQ branches and optional Warburg diffusion.
    FitEcm {
        #[arg(short = 'i', long = "input", required = true, num_args = 1..)]
        input: Vec<PathBuf>,
        #[arg(
            long,
            help = "Circuit: R_QR, R_CR, R_QR_CR, R_CR_CR, or R_QR_QR; append _W for Warburg"
        )]
        model: String,
        #[command(flatten)]
        initial: EcmInitialArgs,
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
        keep_positive_imag: bool,
        #[arg(long)]
        include_correlation_matrix: bool,
        #[command(flatten)]
        batch: BatchArgs,
    },
}

#[derive(Serialize)]
struct DrtSummary {
    lambda: f64,
    tau_min: f64,
    tau_max: f64,
    n_tau: usize,
    tau_grid: TauGridMode,
    basis: DrtBasis,
    shape_control: ShapeControl,
    shape_coefficient: f64,
    epsilon: Option<f64>,
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
    solver: SolverReport,
    constraints: DrtConstraintConfig,
}

#[derive(Serialize)]
struct FitParamsSummary {
    model: String,
    parameters: BTreeMap<String, f64>,
    weighted_sse: f64,
    mean_weighted_chi_square: f64,
    reduced_chi_square: f64,
    relative_rmse_percent: f64,
    rmse_real: f64,
    rmse_imag: f64,
    rmse_magnitude: f64,
    parameter_std_errors: Option<BTreeMap<String, f64>>,
    parameter_rel_std_error_percent: Option<BTreeMap<String, f64>>,
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
        Commands::Clean {
            input,
            lenient,
            imag_sign,
            keep_positive_imag,
            batch,
        } => run_clean(input, lenient, imag_sign, keep_positive_imag, batch),
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
            basis,
            shape_control,
            shape_coefficient,
            regularization_order,
            flip_imag,
            keep_positive_imag,
            nonnegative,
            fit_inductance,
            credible_intervals,
            solver_max_iterations,
            solver_tolerance,
            allow_negative_r_inf,
            nonnegative_inductance,
            compare_matlab_drt,
            compare_matlab_regression,
            batch,
        } => run_drt_batch(
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
            basis,
            shape_control,
            shape_coefficient,
            regularization_order,
            flip_imag,
            keep_positive_imag,
            nonnegative,
            fit_inductance,
            credible_intervals,
            solver_max_iterations,
            solver_tolerance,
            allow_negative_r_inf,
            nonnegative_inductance,
            compare_matlab_drt,
            compare_matlab_regression,
            batch,
        ),
        Commands::FitEcm {
            input,
            model,
            initial,
            auto_init,
            weight,
            max_iter,
            tol,
            flip_imag,
            keep_positive_imag,
            include_correlation_matrix,
            batch,
        } => run_fit_ecm_batch(
            input,
            model,
            initial.into(),
            auto_init,
            weight,
            max_iter,
            tol,
            flip_imag,
            keep_positive_imag,
            include_correlation_matrix,
            batch,
        ),
    }
}

fn run_clean(
    inputs: Vec<PathBuf>,
    lenient: bool,
    imag_sign: ImagSignPolicy,
    keep_positive_imag: bool,
    batch: CleanBatchArgs,
) -> Result<()> {
    let options = BatchOptions {
        jobs: batch.jobs.unwrap_or_else(|| default_jobs(inputs.len())),
        fail_fast: batch.fail_fast,
        overwrite: batch.overwrite,
        resume: false,
        out_root: batch
            .out_root
            .ok_or_else(|| anyhow::anyhow!("output root not specified"))?,
        output_suffix: "cleaned".to_string(),
    };
    let report = run_batch(&inputs, &options, |input, output| {
        clean_file_to(
            input,
            &CleanOptions {
                lenient,
                imag_sign,
                drop_positive_imag: !keep_positive_imag,
                out_root: None,
            },
            output,
        )?;
        Ok(())
    })?;
    finish_batch(report)
}

#[allow(clippy::too_many_arguments)]
fn run_drt_batch(
    inputs: Vec<PathBuf>,
    lambda: f64,
    auto_lambda: bool,
    lambda_min: f64,
    lambda_max: f64,
    n_lambda: usize,
    tau_min: Option<f64>,
    tau_max: Option<f64>,
    n_tau: usize,
    tau_grid: TauGridMode,
    basis: DrtBasis,
    shape_control: ShapeControl,
    shape_coefficient: f64,
    regularization_order: usize,
    flip_imag: bool,
    keep_positive_imag: bool,
    nonnegative: bool,
    fit_inductance: bool,
    credible_intervals: bool,
    solver_max_iterations: usize,
    solver_tolerance: f64,
    allow_negative_r_inf: bool,
    nonnegative_inductance: bool,
    compare_matlab_drt: Option<PathBuf>,
    compare_matlab_regression: Option<PathBuf>,
    batch: BatchArgs,
) -> Result<()> {
    let options = batch_options(&batch, "drt", inputs.len());
    let configuration = json!({
        "input_policy": {"flip_imag": flip_imag, "drop_positive_imag": !keep_positive_imag},
        "lambda": lambda, "auto_lambda": auto_lambda, "lambda_min": lambda_min,
        "lambda_max": lambda_max, "n_lambda": n_lambda,
        "tau_min": tau_min, "tau_max": tau_max, "n_tau": n_tau, "tau_grid": tau_grid,
        "basis": basis, "shape_control": shape_control, "shape_coefficient": shape_coefficient,
        "regularization_order": regularization_order, "nonnegative": nonnegative,
        "fit_inductance": fit_inductance, "credible_intervals": credible_intervals,
        "solver_max_iterations": solver_max_iterations, "solver_tolerance": solver_tolerance,
        "constraints": {"gamma_nonnegative": true, "r_inf_nonnegative": !allow_negative_r_inf,
            "inductance_nonnegative": nonnegative_inductance},
        "compare_matlab_drt": compare_matlab_drt, "compare_matlab_regression": compare_matlab_regression,
    });
    let report = run_batch_with_resume(
        &inputs,
        &options,
        |input, output| {
            let expected = RunManifest::new("drt", input, configuration.clone())?;
            verify_resume(output, &expected)
        },
        |input, output| {
            let manifest = RunManifest::new("drt", input, configuration.clone())?;
            execute_manifested(output, manifest, || {
                run_drt(
                    input.to_path_buf(),
                    lambda,
                    auto_lambda,
                    lambda_min,
                    lambda_max,
                    n_lambda,
                    tau_min,
                    tau_max,
                    n_tau,
                    tau_grid,
                    basis,
                    shape_control,
                    shape_coefficient,
                    regularization_order,
                    flip_imag,
                    keep_positive_imag,
                    nonnegative,
                    fit_inductance,
                    credible_intervals,
                    solver_max_iterations,
                    solver_tolerance,
                    allow_negative_r_inf,
                    nonnegative_inductance,
                    compare_matlab_drt.clone(),
                    compare_matlab_regression.clone(),
                    Some(output.to_path_buf()),
                )?;
                collect_output_files(output)
            })
        },
    )?;
    finish_batch(report)
}

#[allow(clippy::too_many_arguments)]
fn run_fit_ecm_batch(
    inputs: Vec<PathBuf>,
    model_name: String,
    initial: PartialEcmInit,
    auto_init: bool,
    weight: Weighting,
    max_iter: usize,
    tol: f64,
    flip_imag: bool,
    keep_positive_imag: bool,
    include_correlation_matrix: bool,
    batch: BatchArgs,
) -> Result<()> {
    let model: EcmModelSpec = model_name.parse()?;
    let options = batch_options(&batch, "fit_ecm", inputs.len());
    let configuration = json!({
        "input_policy": {"flip_imag": flip_imag, "drop_positive_imag": !keep_positive_imag},
        "model": model.canonical_name(),
        "initial": {
            "rs": initial.rs,
            "r1": initial.r1, "c1": initial.c1, "q1": initial.q1, "n1": initial.n1,
            "r2": initial.r2, "c2": initial.c2, "q2": initial.q2, "n2": initial.n2,
            "warburg_sigma": initial.warburg_sigma,
        },
        "auto_init": auto_init, "weight": weight, "max_iterations": max_iter,
        "tolerance": tol, "include_correlation_matrix": include_correlation_matrix,
    });
    let report = run_batch_with_resume(
        &inputs,
        &options,
        |input, output| {
            let expected = RunManifest::new("fit-ecm", input, configuration.clone())?;
            verify_resume(output, &expected)
        },
        |input, output| {
            let manifest = RunManifest::new("fit-ecm", input, configuration.clone())?;
            execute_manifested(output, manifest, || {
                run_fit_ecm(
                    input.to_path_buf(),
                    model.clone(),
                    initial.clone(),
                    auto_init,
                    weight,
                    max_iter,
                    tol,
                    flip_imag,
                    keep_positive_imag,
                    include_correlation_matrix,
                    Some(output.to_path_buf()),
                )?;
                collect_output_files(output)
            })
        },
    )?;
    finish_batch(report)
}

fn batch_options(args: &BatchArgs, process_type: &str, input_count: usize) -> BatchOptions {
    BatchOptions {
        jobs: args.jobs.unwrap_or_else(|| default_jobs(input_count)),
        fail_fast: args.fail_fast,
        overwrite: args.overwrite,
        resume: args.resume,
        out_root: args.out_root.clone().unwrap_or_else(|| PathBuf::from("")),
        output_suffix: process_type.to_string(),
    }
}

fn finish_batch(report: BatchReport) -> Result<()> {
    for item in &report.items {
        match item.status {
            BatchStatus::Success => println!(
                "ok      {} -> {}",
                item.input_path.display(),
                item.output_dir.display()
            ),
            BatchStatus::Resumed => println!(
                "resumed {} -> {}",
                item.input_path.display(),
                item.output_dir.display()
            ),
            BatchStatus::NotProcessed => println!(
                "not processed {} -> {}",
                item.input_path.display(),
                item.output_dir.display()
            ),
            BatchStatus::Failed => eprintln!(
                "failed  {}: {}",
                item.input_path.display(),
                item.error.as_deref().unwrap_or("unknown error")
            ),
        }
    }
    println!(
        "batch: {} succeeded, {} failed, {} resumed, {} not processed",
        report.succeeded, report.failed, report.resumed, report.not_processed
    );
    if report.failed > 0 {
        bail!("batch completed with {} failed input(s)", report.failed);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
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
    basis: DrtBasis,
    shape_control: ShapeControl,
    shape_coefficient: f64,
    regularization_order: usize,
    flip_imag: bool,
    keep_positive_imag: bool,
    nonnegative: bool,
    fit_inductance: bool,
    credible_intervals: bool,
    solver_max_iterations: usize,
    solver_tolerance: f64,
    allow_negative_r_inf: bool,
    nonnegative_inductance: bool,
    compare_matlab_drt: Option<PathBuf>,
    compare_matlab_regression: Option<PathBuf>,
    out: Option<PathBuf>,
) -> Result<()> {
    let data = read_analysis_data(&input, flip_imag, keep_positive_imag, "DRT")?;
    let out = out.unwrap_or_else(|| default_output_dir(&input, "drt"));
    fs::create_dir_all(&out).with_context(|| format!("failed to create {}", out.display()))?;

    let mut settings = DrtSettings {
        lambda,
        tau_min,
        tau_max,
        n_tau,
        tau_grid,
        basis,
        shape_control,
        shape_coefficient,
        fit_inductance,
        regularization_order,
        nonnegative,
        credible_intervals,
        solver: DrtSolverOptions {
            max_iterations: solver_max_iterations,
            tolerance: solver_tolerance,
            constraints: DrtConstraintConfig {
                gamma_nonnegative: true,
                r_inf_nonnegative: !allow_negative_r_inf,
                inductance_nonnegative: nonnegative_inductance,
            },
        },
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
        basis: result.settings_used.basis,
        shape_control: result.settings_used.shape_control,
        shape_coefficient: result.settings_used.shape_coefficient,
        epsilon: result.settings_used.epsilon,
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
        solver: result.solver_report.clone(),
        constraints: result.settings_used.constraints,
    };
    fs::write(
        out.join("residual_summary.json"),
        serde_json::to_string_pretty(&summary)?,
    )?;
    fs::write(
        out.join("solver_report.json"),
        serde_json::to_string_pretty(&result.solver_report)?,
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
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_fit_ecm(
    input: PathBuf,
    model: EcmModelSpec,
    partial_initial: PartialEcmInit,
    auto_init: bool,
    weight: Weighting,
    max_iter: usize,
    tol: f64,
    flip_imag: bool,
    keep_positive_imag: bool,
    include_correlation_matrix: bool,
    out: Option<PathBuf>,
) -> Result<()> {
    let data = read_analysis_data(&input, flip_imag, keep_positive_imag, "ECM fitting")?;
    let out = out.unwrap_or_else(|| default_output_dir(&input, "ecm"));
    fs::create_dir_all(&out).with_context(|| format!("failed to create {}", out.display()))?;

    let initial: EcmParams = complete_initial_ecm(&data, &model, partial_initial, auto_init)?;
    let result = fit_ecm(
        &data,
        &EcmFitSettings {
            model: model.clone(),
            initial,
            weight,
            max_iter,
            tol,
        },
    )?;

    write_impedance_csv(&out.join("fitted_impedance.csv"), &data, &result.z_fit)?;
    let model_name = model.canonical_name();
    let plot_title = format!("{} Nyquist Fit", model_name);
    write_nyquist_svg(
        &out.join("nyquist_fit.svg"),
        &plot_title,
        &data.z_real,
        &data.z_imag,
        &result.z_fit,
    )?;
    let labels = model.parameter_labels();
    let summary = FitParamsSummary {
        model: model_name,
        parameters: labeled_parameters(&labels, &result.params.values()),
        weighted_sse: result.weighted_sse,
        mean_weighted_chi_square: result.mean_weighted_chi_square,
        reduced_chi_square: result.reduced_chi_square,
        relative_rmse_percent: result.metrics.relative_rmse_percent,
        rmse_real: result.metrics.rmse_real,
        rmse_imag: result.metrics.rmse_imag,
        rmse_magnitude: result.metrics.rmse_magnitude,
        parameter_std_errors: result
            .parameter_std_errors
            .as_deref()
            .map(|values| labeled_parameters(&labels, values)),
        parameter_rel_std_error_percent: result
            .parameter_rel_std_error_percent
            .as_deref()
            .map(|values| labeled_parameters(&labels, values)),
        parameter_correlation_labels: labels,
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
    Ok(())
}

fn read_analysis_data(
    input: &Path,
    flip_imag: bool,
    keep_positive_imag: bool,
    analysis: &str,
) -> Result<EisData> {
    let mut data = read_eis_with_cleaning(input, false)?;
    if flip_imag {
        data.flip_imag();
    }
    if keep_positive_imag {
        data.warn_if_imag_mostly_positive();
    } else {
        let original_count = data.len();
        let removed = data.drop_positive_imag()?;
        println!(
            "positive-imaginary filter ({analysis}): removed {removed} of {original_count} point(s) from {}; use --keep-positive-imag to disable this filter",
            input.display()
        );
    }
    Ok(data)
}

fn labeled_parameters(labels: &[String], values: &[f64]) -> BTreeMap<String, f64> {
    labels.iter().cloned().zip(values.iter().copied()).collect()
}

fn write_drt_svgs(
    out: &Path,
    data: &electrochem_tools::types::EisData,
    result: &electrochem_tools::drt::DrtResult,
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

fn format_matlab_exp(value: f64) -> String {
    let raw = format!("{value:.6e}");
    let Some((mantissa, exponent)) = raw.split_once('e') else {
        return raw;
    };
    let exp_value = exponent.parse::<i32>().unwrap_or(0);
    format!("{mantissa}e{exp_value:+03}")
}
