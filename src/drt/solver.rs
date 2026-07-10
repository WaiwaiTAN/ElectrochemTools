use crate::regularization::{
    solve_tikhonov_active_set_with_penalty, solve_tikhonov_general_with_penalty,
};
use anyhow::Result;
use nalgebra::{DMatrix, DVector};

pub fn solve_coefficients(
    matrix: &DMatrix<f64>,
    observations: &DVector<f64>,
    lambda: f64,
    penalty: &DMatrix<f64>,
    nonnegative: bool,
) -> Result<DVector<f64>> {
    if nonnegative {
        let lower_bound_zero = vec![true; matrix.ncols()];
        solve_tikhonov_active_set_with_penalty(
            matrix,
            observations,
            lambda,
            penalty,
            &lower_bound_zero,
            1_000,
            1.0e-9,
        )
    } else {
        solve_tikhonov_general_with_penalty(matrix, observations, lambda, penalty)
    }
}
