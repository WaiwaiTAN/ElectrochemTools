use crate::regularization::{
    SolverReport, solve_tikhonov_active_set_with_penalty_report,
    solve_tikhonov_general_with_penalty_report,
};
use anyhow::{Result, bail};
use nalgebra::{DMatrix, DVector};
use serde::Serialize;

#[derive(Debug, Clone, Copy, Serialize)]
pub struct DrtConstraintConfig {
    pub gamma_nonnegative: bool,
    pub r_inf_nonnegative: bool,
    pub inductance_nonnegative: bool,
}

impl Default for DrtConstraintConfig {
    fn default() -> Self {
        Self {
            gamma_nonnegative: true,
            r_inf_nonnegative: true,
            inductance_nonnegative: false,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct DrtSolverOptions {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub constraints: DrtConstraintConfig,
}

impl Default for DrtSolverOptions {
    fn default() -> Self {
        Self {
            max_iterations: 1_000,
            tolerance: 1.0e-9,
            constraints: DrtConstraintConfig::default(),
        }
    }
}

pub struct CoefficientSolution {
    pub coefficients: DVector<f64>,
    pub report: SolverReport,
}

pub fn solve_coefficients(
    matrix: &DMatrix<f64>,
    observations: &DVector<f64>,
    lambda: f64,
    penalty: &DMatrix<f64>,
    nonnegative: bool,
    fit_inductance: bool,
    options: DrtSolverOptions,
) -> Result<CoefficientSolution> {
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        bail!("solver tolerance must be finite and strictly positive");
    }
    if nonnegative {
        let mut lower_bound_zero = vec![options.constraints.gamma_nonnegative; matrix.ncols()];
        if fit_inductance {
            lower_bound_zero[0] = options.constraints.inductance_nonnegative;
            lower_bound_zero[1] = options.constraints.r_inf_nonnegative;
        } else {
            lower_bound_zero[0] = options.constraints.r_inf_nonnegative;
        }
        let (coefficients, report) = solve_tikhonov_active_set_with_penalty_report(
            matrix,
            observations,
            lambda,
            penalty,
            &lower_bound_zero,
            options.max_iterations,
            options.tolerance,
        )?;
        if !report.converged {
            bail!(
                "DRT active-set solver did not converge after {} iterations (KKT violation {:.3e})",
                report.iterations,
                report.kkt_violation
            );
        }
        Ok(CoefficientSolution {
            coefficients,
            report,
        })
    } else {
        let (coefficients, report) =
            solve_tikhonov_general_with_penalty_report(matrix, observations, lambda, penalty)?;
        if !report.converged {
            bail!("unconstrained DRT solve failed the gradient convergence check");
        }
        Ok(CoefficientSolution {
            coefficients,
            report,
        })
    }
}
