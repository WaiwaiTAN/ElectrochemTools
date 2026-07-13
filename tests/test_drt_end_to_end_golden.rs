use electrochem_tools::drt::{
    DrtConstraintConfig, DrtSettings, DrtSolverOptions, TauGridMode, solve_drt,
};
use electrochem_tools::eis_io::read_eis_with_cleaning;
use serde::Deserialize;
use std::path::{Path, PathBuf};

#[derive(Deserialize)]
struct Summary {
    lambda: f64,
    regularization_order: usize,
    fit_inductance: bool,
    objective_value: f64,
    #[serde(rename = "R_inf")]
    r_inf: f64,
    inductance: f64,
    polarization_resistance: f64,
    constraints: Constraints,
}

#[derive(Deserialize)]
struct Constraints {
    gamma_nonnegative: bool,
    r_inf_nonnegative: bool,
    inductance_mode: String,
}

fn numbers(path: &Path) -> Vec<f64> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .unwrap();
    reader
        .records()
        .flat_map(|row| {
            row.unwrap()
                .iter()
                .map(|value| value.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect()
}

fn columns(path: &Path) -> Vec<[f64; 2]> {
    numbers(path)
        .chunks_exact(2)
        .map(|row| [row[0], row[1]])
        .collect()
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
        .map(|(a, e)| (a - e).abs())
        .fold(0.0, f64::max);
    let difference_l2 = actual
        .iter()
        .zip(expected)
        .map(|(a, e)| (a - e).powi(2))
        .sum::<f64>()
        .sqrt();
    let expected_l2 = expected
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    (max_absolute, difference_l2 / expected_l2.max(1.0e-300))
}

fn compare_case(name: &str) {
    let root = PathBuf::from("tests/golden/drttools/end_to_end").join(name);
    let summary: Summary =
        serde_json::from_slice(&std::fs::read(root.join("summary.json")).unwrap()).unwrap();
    assert!(summary.constraints.gamma_nonnegative);
    assert!(summary.constraints.r_inf_nonnegative);
    assert_eq!(
        summary.constraints.inductance_mode,
        if summary.fit_inductance {
            "nonnegative"
        } else {
            "fixed_zero"
        }
    );

    let data = read_eis_with_cleaning(Path::new("examples/data/eis_cleaned.csv"), false).unwrap();
    let result = solve_drt(
        &data,
        &DrtSettings {
            lambda: summary.lambda,
            tau_min: None,
            tau_max: None,
            n_tau: data.len(),
            tau_grid: TauGridMode::Drttools,
            basis: Default::default(),
            shape_control: Default::default(),
            shape_coefficient: 0.5,
            fit_inductance: summary.fit_inductance,
            regularization_order: summary.regularization_order,
            nonnegative: true,
            credible_intervals: false,
            solver: DrtSolverOptions {
                max_iterations: 1_000,
                tolerance: 1.0e-9,
                constraints: DrtConstraintConfig {
                    gamma_nonnegative: true,
                    r_inf_nonnegative: true,
                    inductance_nonnegative: summary.fit_inductance,
                },
            },
        },
    )
    .unwrap();

    let frequency = numbers(&root.join("frequency.csv"));
    let tau = numbers(&root.join("tau.csv"));
    assert_eq!(frequency, data.frequency_hz);
    for (index, (&actual, &expected)) in result.tau.iter().zip(&tau).enumerate() {
        mixed_close(actual, expected, 1e-12, 1e-12, &format!("tau[{index}]"));
    }
    mixed_close(result.r_inf, summary.r_inf, 1e-11, 1e-9, "R_inf");
    mixed_close(
        result.inductance,
        summary.inductance,
        1e-12,
        1e-9,
        "inductance",
    );
    mixed_close(
        result.polarization_resistance,
        summary.polarization_resistance,
        1e-10,
        1e-9,
        "polarization resistance",
    );
    mixed_close(
        2.0 * result.solver_report.objective,
        summary.objective_value,
        1e-10,
        1e-9,
        "objective",
    );

    let expected_coefficients = numbers(&root.join("coefficients.csv"));
    let mut coefficients = vec![result.inductance, result.r_inf];
    coefficients.extend_from_slice(&result.gamma);
    let (coefficient_max, coefficient_l2) = vector_errors(&coefficients, &expected_coefficients);
    let expected_gamma = numbers(&root.join("gamma.csv"));
    let (gamma_max, gamma_l2) = vector_errors(&result.gamma, &expected_gamma);
    assert!(
        coefficient_max <= 1e-9 && coefficient_l2 <= 1e-10,
        "coefficients: max_abs={coefficient_max:.3e}, l2_rel={coefficient_l2:.3e}"
    );
    assert!(
        gamma_max <= 1e-9 && gamma_l2 <= 1e-10,
        "gamma: max_abs={gamma_max:.3e}, l2_rel={gamma_l2:.3e}"
    );

    let expected_z = columns(&root.join("reconstructed_impedance.csv"));
    let actual_z = result
        .z_fit
        .iter()
        .flat_map(|z| [z.re, z.im])
        .collect::<Vec<_>>();
    let expected_z_flat = expected_z.iter().flat_map(|z| *z).collect::<Vec<_>>();
    let (z_max, z_l2) = vector_errors(&actual_z, &expected_z_flat);
    assert!(
        z_max <= 1e-10 && z_l2 <= 1e-11,
        "reconstructed impedance: max_abs={z_max:.3e}, l2_rel={z_l2:.3e}"
    );
    println!(
        "{name}: coeff max={coefficient_max:.3e} l2={coefficient_l2:.3e}; gamma max={gamma_max:.3e} l2={gamma_l2:.3e}; Z max={z_max:.3e} l2={z_l2:.3e}; objective rel={:.3e}",
        ((2.0 * result.solver_report.objective - summary.objective_value)
            / summary.objective_value)
            .abs()
    );
}

#[test]
fn matlab_drttools_end_to_end_cases_match() {
    for name in [
        "real_order1_no_inductance",
        "real_order1_with_inductance",
        "real_order2_no_inductance",
    ] {
        compare_case(name);
    }
}
