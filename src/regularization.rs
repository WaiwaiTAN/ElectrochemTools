use anyhow::{Context, Result, bail};
use nalgebra::{DMatrix, DVector};

pub fn solve_tikhonov(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    n_gamma: usize,
    order: usize,
) -> Result<DVector<f64>> {
    if a.nrows() != b.len() {
        bail!("A row count must match b length");
    }
    if a.ncols() != n_gamma + 1 {
        bail!("A column count must be 1 + number of gamma values");
    }
    if lambda < 0.0 || !lambda.is_finite() {
        bail!("lambda must be finite and non-negative");
    }

    let ata = a.transpose() * a;
    let atb = a.transpose() * b;
    let penalty = penalty_matrix(n_gamma, order)?;
    let system = ata + penalty.scale(lambda);

    if let Some(cholesky) = system.clone().cholesky() {
        return Ok(cholesky.solve(&atb));
    }
    system
        .lu()
        .solve(&atb)
        .context("Tikhonov system is singular; try a larger lambda")
}

pub fn penalty_matrix(n_gamma: usize, order: usize) -> Result<DMatrix<f64>> {
    let n_params = n_gamma + 1;
    let mut l = match order {
        0 => DMatrix::<f64>::zeros(n_gamma, n_params),
        1 => {
            if n_gamma < 2 {
                bail!("first-order regularization requires at least two gamma values");
            }
            DMatrix::<f64>::zeros(n_gamma - 1, n_params)
        }
        2 => {
            if n_gamma < 3 {
                bail!("second-order regularization requires at least three gamma values");
            }
            DMatrix::<f64>::zeros(n_gamma - 2, n_params)
        }
        _ => bail!("regularization order must be 0, 1, or 2"),
    };

    match order {
        0 => {
            for row in 0..n_gamma {
                l[(row, row + 1)] = 1.0;
            }
        }
        1 => {
            for row in 0..(n_gamma - 1) {
                l[(row, row + 1)] = -1.0;
                l[(row, row + 2)] = 1.0;
            }
        }
        2 => {
            for row in 0..(n_gamma - 2) {
                l[(row, row + 1)] = 1.0;
                l[(row, row + 2)] = -2.0;
                l[(row, row + 3)] = 1.0;
            }
        }
        _ => unreachable!(),
    }
    Ok(l.transpose() * l)
}
