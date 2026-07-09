use ElectrochemTools::drt::{DrtSettings, solve_drt};
use ElectrochemTools::ecm::RqrParams;
use ElectrochemTools::eis_io::{read_eis, write_impedance_csv};
use ElectrochemTools::fit::{
    PartialRqrInit, RqrFitSettings, Weighting, complete_initial_params, fit_rqr,
};
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
        tau_min: Option<f64>,
        #[arg(long)]
        tau_max: Option<f64>,
        #[arg(long, default_value_t = 100)]
        n_tau: usize,
        #[arg(long, default_value_t = 1)]
        regularization_order: usize,
        #[arg(long)]
        flip_imag: bool,
        #[arg(long)]
        out: PathBuf,
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
        out: PathBuf,
    },
}

#[derive(Serialize)]
struct DrtSummary {
    lambda: f64,
    tau_min: f64,
    tau_max: f64,
    n_tau: usize,
    regularization_order: usize,
    n_points: usize,
    r_inf: f64,
    rmse_real: f64,
    rmse_imag: f64,
    rmse_magnitude: f64,
    relative_rmse_percent: f64,
    note: String,
}

#[derive(Serialize)]
struct FitParamsSummary {
    model: String,
    rs: f64,
    rct: f64,
    q: f64,
    n: f64,
    relative_rmse_percent: f64,
    rmse_real: f64,
    rmse_imag: f64,
    rmse_magnitude: f64,
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
            tau_min,
            tau_max,
            n_tau,
            regularization_order,
            flip_imag,
            out,
        } => run_drt(
            input,
            lambda,
            tau_min,
            tau_max,
            n_tau,
            regularization_order,
            flip_imag,
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
            out,
        } => run_fit_ecm(
            input, model, rs, rct, q, n, auto_init, weight, max_iter, tol, flip_imag, out,
        ),
    }
}

fn run_drt(
    input: PathBuf,
    lambda: f64,
    tau_min: Option<f64>,
    tau_max: Option<f64>,
    n_tau: usize,
    regularization_order: usize,
    flip_imag: bool,
    out: PathBuf,
) -> Result<()> {
    let mut data = read_eis(&input)?;
    if flip_imag {
        data.flip_imag();
    }
    data.warn_if_imag_mostly_positive();
    fs::create_dir_all(&out).with_context(|| format!("failed to create {}", out.display()))?;

    let result = solve_drt(
        &data,
        &DrtSettings {
            lambda,
            tau_min,
            tau_max,
            n_tau,
            regularization_order,
        },
    )?;

    let mut gamma_writer = csv::Writer::from_path(out.join("gamma.csv"))?;
    gamma_writer.write_record(["tau", "log10_tau", "gamma"])?;
    for (&tau, &gamma) in result.tau.iter().zip(&result.gamma) {
        gamma_writer.serialize((tau, tau.log10(), gamma))?;
    }
    gamma_writer.flush()?;

    write_impedance_csv(
        &out.join("reconstructed_impedance.csv"),
        &data,
        &result.z_fit,
    )?;
    let summary = DrtSummary {
        lambda: result.settings_used.lambda,
        tau_min: result.settings_used.tau_min,
        tau_max: result.settings_used.tau_max,
        n_tau: result.settings_used.n_tau,
        regularization_order: result.settings_used.regularization_order,
        n_points: data.len(),
        r_inf: result.r_inf,
        rmse_real: result.metrics.rmse_real,
        rmse_imag: result.metrics.rmse_imag,
        rmse_magnitude: result.metrics.rmse_magnitude,
        relative_rmse_percent: result.metrics.relative_rmse_percent,
        note: "unconstrained Tikhonov DRT MVP; gamma nonnegativity not enforced".to_string(),
    };
    fs::write(
        out.join("residual_summary.json"),
        serde_json::to_string_pretty(&summary)?,
    )?;
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
    out: PathBuf,
) -> Result<()> {
    if model.to_ascii_uppercase() != "R_QR" {
        bail!("unsupported model '{}'; currently supported: R_QR", model);
    }

    let mut data = read_eis(&input)?;
    if flip_imag {
        data.flip_imag();
    }
    data.warn_if_imag_mostly_positive();
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
    let summary = FitParamsSummary {
        model: "R_QR".to_string(),
        rs: result.params.rs,
        rct: result.params.rct,
        q: result.params.q,
        n: result.params.n,
        relative_rmse_percent: result.metrics.relative_rmse_percent,
        rmse_real: result.metrics.rmse_real,
        rmse_imag: result.metrics.rmse_imag,
        rmse_magnitude: result.metrics.rmse_magnitude,
        weight: result.weight,
        n_points: data.len(),
        iterations: result.iterations,
        converged: result.converged,
    };
    fs::write(
        out.join("fit_params.json"),
        serde_json::to_string_pretty(&summary)?,
    )?;
    println!("ECM fitting complete: {}", out.display());
    Ok(())
}
