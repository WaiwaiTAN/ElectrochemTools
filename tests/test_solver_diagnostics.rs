use electrochem_tools::drt::solver::{DrtConstraintConfig, DrtSolverOptions, solve_coefficients};
use electrochem_tools::regularization::solve_tikhonov_active_set_with_penalty_report;
use nalgebra::{DMatrix, DVector};

#[test]
fn constraints_are_configured_by_physical_parameter_role() {
    let matrix = DMatrix::<f64>::identity(3, 3);
    let observations = DVector::from_vec(vec![-2.0, -1.0, 3.0]);
    let penalty = DMatrix::<f64>::zeros(3, 3);
    let default = solve_coefficients(
        &matrix,
        &observations,
        0.0,
        &penalty,
        true,
        true,
        DrtSolverOptions::default(),
    )
    .unwrap();
    assert!((default.coefficients[0] + 2.0).abs() < 1.0e-12); // inductance is free
    assert_eq!(default.coefficients[1], 0.0); // R_inf is bounded
    assert!((default.coefficients[2] - 3.0).abs() < 1.0e-12); // gamma is bounded
    assert!(default.report.converged);
    assert_eq!(default.report.active_constraints, 1);
    assert!(default.report.condition_estimate.is_some());

    let alternate = solve_coefficients(
        &matrix,
        &observations,
        0.0,
        &penalty,
        true,
        true,
        DrtSolverOptions {
            constraints: DrtConstraintConfig {
                gamma_nonnegative: true,
                r_inf_nonnegative: false,
                inductance_nonnegative: true,
            },
            ..DrtSolverOptions::default()
        },
    )
    .unwrap();
    assert_eq!(alternate.coefficients[0], 0.0);
    assert!((alternate.coefficients[1] + 1.0).abs() < 1.0e-12);
}

#[test]
fn maximum_iteration_exit_is_structured_and_not_converged() {
    let matrix = DMatrix::<f64>::identity(2, 2);
    let observations = DVector::from_vec(vec![1.0, -1.0]);
    let penalty = DMatrix::<f64>::zeros(2, 2);
    let (_, report) = solve_tikhonov_active_set_with_penalty_report(
        &matrix,
        &observations,
        0.0,
        &penalty,
        &[true, true],
        0,
        1.0e-12,
    )
    .unwrap();
    assert!(!report.converged);
    assert_eq!(report.iterations, 0);
    assert!(
        report
            .warning
            .as_deref()
            .unwrap()
            .contains("maximum iterations")
    );

    let error = solve_coefficients(
        &matrix,
        &observations,
        0.0,
        &penalty,
        true,
        false,
        DrtSolverOptions {
            max_iterations: 0,
            ..DrtSolverOptions::default()
        },
    )
    .err()
    .unwrap();
    assert!(error.to_string().contains("did not converge"));
}
