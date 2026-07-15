use electrochem_tools::drt::discretization::{DrtDiscretization, GaussianDiscretization};
use electrochem_tools::drt::kernel::assemble_combined_from_kernels;
use electrochem_tools::drt::solver::solve_coefficients;
use electrochem_tools::drt::{
    DrtBasis, DrtConstraintConfig, DrtSettings, DrtSolverOptions, ShapeControl, TauGridMode,
    make_drttools_tau_grid, solve_drt,
};
use electrochem_tools::eis_io::read_eis_with_cleaning;
use nalgebra::DMatrix;
use serde::Deserialize;
use std::path::{Path, PathBuf};

#[derive(Deserialize)]
struct Summary {
    lambda: f64,
    regularization_order: usize,
    shape_coefficient: f64,
    epsilon: f64,
    objective_value: f64,
    #[serde(rename = "R_inf")]
    r_inf: f64,
    polarization_resistance: f64,
}

fn root() -> PathBuf {
    PathBuf::from("tests/golden/drttools/end_to_end/gaussian_simple_run")
}

fn matrix(path: impl AsRef<Path>) -> DMatrix<f64> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .unwrap();
    let rows = reader
        .records()
        .map(|record| {
            record
                .unwrap()
                .iter()
                .map(|value| value.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let columns = rows.first().unwrap().len();
    DMatrix::from_row_iterator(rows.len(), columns, rows.into_iter().flatten())
}

fn values(name: &str) -> Vec<f64> {
    matrix(root().join(name)).iter().copied().collect()
}

fn mixed_close(actual: f64, expected: f64, absolute: f64, relative: f64, label: &str) {
    let tolerance = absolute.max(relative * expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{label}: actual={actual:.16e}, expected={expected:.16e}, error={:.3e}, tolerance={tolerance:.3e}",
        (actual - expected).abs()
    );
}

fn vector_errors(actual: &[f64], expected: &[f64]) -> (f64, f64) {
    assert_eq!(actual.len(), expected.len());
    let max_absolute = actual
        .iter()
        .zip(expected)
        .map(|(actual, expected)| (actual - expected).abs())
        .fold(0.0, f64::max);
    let difference_l2 = actual
        .iter()
        .zip(expected)
        .map(|(actual, expected)| (actual - expected).powi(2))
        .sum::<f64>()
        .sqrt();
    let expected_l2 = expected
        .iter()
        .map(|expected| expected * expected)
        .sum::<f64>()
        .sqrt();
    (max_absolute, difference_l2 / expected_l2.max(1.0e-300))
}

#[test]
fn gaussian_simple_run_matches_matlab_drttools() {
    let summary: Summary =
        serde_json::from_slice(&std::fs::read(root().join("summary.json")).unwrap()).unwrap();
    assert_eq!(summary.regularization_order, 1);
    let data = read_eis_with_cleaning(Path::new("tests/fixtures/eis_cleaned.csv"), false).unwrap();
    let tau = make_drttools_tau_grid(&data).unwrap();
    let discretization =
        GaussianDiscretization::new(tau.clone(), ShapeControl::Fwhm, summary.shape_coefficient)
            .unwrap();
    mixed_close(
        discretization.epsilon().unwrap(),
        summary.epsilon,
        1.0e-13,
        1.0e-13,
        "epsilon",
    );

    let kernels = discretization.kernel_matrices(&data.frequency_hz).unwrap();
    let matlab_a_re = matrix(root().join("A_re.csv"));
    let matlab_a_im = matrix(root().join("A_im.csv"));
    for row in 0..data.len() {
        mixed_close(matlab_a_re[(row, 1)], 1.0, 0.0, 0.0, "R column");
        for column in 0..tau.len() {
            mixed_close(
                kernels.real[(row, column)],
                matlab_a_re[(row, column + 2)],
                2.0e-10,
                2.0e-9,
                "A_re",
            );
            mixed_close(
                kernels.imaginary[(row, column)],
                matlab_a_im[(row, column + 2)],
                2.0e-10,
                2.0e-9,
                "A_im",
            );
        }
    }
    let penalty = discretization.penalty(1, 2).unwrap();
    let matlab_penalty = matrix(root().join("M_1.csv"));
    for (index, (&actual, &expected)) in penalty.iter().zip(matlab_penalty.iter()).enumerate() {
        mixed_close(actual, expected, 1.0e-12, 1.0e-11, &format!("M_1[{index}]"));
    }

    let (combined, observations) = assemble_combined_from_kernels(&data, &kernels, false);
    let compact_penalty = discretization.penalty(1, 1).unwrap();
    let coefficient_solution = solve_coefficients(
        &combined,
        &observations,
        summary.lambda,
        &compact_penalty,
        true,
        false,
        DrtSolverOptions {
            max_iterations: 1_000,
            tolerance: 1.0e-9,
            constraints: DrtConstraintConfig::default(),
        },
    )
    .unwrap();
    let expected_coefficients = values("coefficients.csv");
    let mut actual_coefficients = vec![0.0];
    actual_coefficients.extend(coefficient_solution.coefficients.iter().copied());
    let (coefficient_max, coefficient_l2) =
        vector_errors(&actual_coefficients, &expected_coefficients);
    assert!(
        coefficient_max <= 2.0e-9 && coefficient_l2 <= 2.0e-9,
        "coefficients: max_abs={coefficient_max:.3e}, l2_rel={coefficient_l2:.3e}"
    );

    let result = solve_drt(
        &data,
        &DrtSettings {
            lambda: summary.lambda,
            tau_min: None,
            tau_max: None,
            n_tau: 100,
            tau_grid: TauGridMode::Logspace,
            basis: DrtBasis::Gaussian,
            shape_control: ShapeControl::Fwhm,
            shape_coefficient: summary.shape_coefficient,
            fit_inductance: false,
            regularization_order: 1,
            nonnegative: true,
            credible_intervals: false,
            bayesian: None,
            solver: DrtSolverOptions::default(),
        },
    )
    .unwrap();
    assert!(matches!(
        result.settings_used.tau_grid,
        TauGridMode::Drttools
    ));
    assert_eq!(result.plot_tau.len(), 10 * result.tau.len());
    mixed_close(
        result.plot_tau[0].log10(),
        result.tau[0].log10() - 0.5,
        1.0e-14,
        1.0e-14,
        "dense plot tau minimum",
    );
    mixed_close(
        result.plot_tau[result.plot_tau.len() - 1].log10(),
        result.tau[result.tau.len() - 1].log10() + 0.5,
        1.0e-14,
        1.0e-14,
        "dense plot tau maximum",
    );
    mixed_close(result.r_inf, summary.r_inf, 2.0e-10, 2.0e-10, "R_inf");
    mixed_close(
        result.polarization_resistance,
        summary.polarization_resistance,
        2.0e-10,
        2.0e-9,
        "polarization resistance",
    );
    mixed_close(
        2.0 * result.solver_report.objective,
        summary.objective_value,
        2.0e-9,
        2.0e-10,
        "objective",
    );

    let expected_gamma = values("gamma.csv");
    let (gamma_max, gamma_l2) = vector_errors(&result.gamma, &expected_gamma);
    assert!(
        gamma_max <= 2.0e-9 && gamma_l2 <= 2.0e-9,
        "mapped gamma: max_abs={gamma_max:.3e}, l2_rel={gamma_l2:.3e}"
    );
    let expected_z = matrix(root().join("reconstructed_impedance.csv"));
    let actual_z = result
        .z_fit
        .iter()
        .flat_map(|value| [value.re, value.im])
        .collect::<Vec<_>>();
    let expected_z = expected_z
        .row_iter()
        .flat_map(|row| [row[0], row[1]])
        .collect::<Vec<_>>();
    let (z_max, z_l2) = vector_errors(&actual_z, &expected_z);
    assert!(
        z_max <= 2.0e-9 && z_l2 <= 2.0e-9,
        "reconstructed impedance: max_abs={z_max:.3e}, l2_rel={z_l2:.3e}"
    );
    println!(
        "gaussian_simple_run: coeff max={coefficient_max:.3e} l2={coefficient_l2:.3e}; gamma max={gamma_max:.3e} l2={gamma_l2:.3e}; Z max={z_max:.3e} l2={z_l2:.3e}; objective rel={:.3e}",
        ((2.0 * result.solver_report.objective - summary.objective_value)
            / summary.objective_value)
            .abs()
    );
}
