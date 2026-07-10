use anyhow::{Result, bail};
use nalgebra::DMatrix;

pub fn piecewise_linear_penalty(
    tau: &[f64],
    order: usize,
    n_unpenalized: usize,
) -> Result<DMatrix<f64>> {
    let n_gamma = tau.len();
    let n_params = n_gamma + n_unpenalized;
    let rows = match order {
        0 => n_gamma,
        1 if n_gamma >= 2 => n_gamma - 1,
        1 => bail!("first-order regularization requires at least two gamma values"),
        2 if n_gamma >= 3 => n_gamma - 2,
        2 => bail!("second-order regularization requires at least three gamma values"),
        _ => bail!("regularization order must be 0, 1, or 2"),
    };
    let mut derivative = DMatrix::<f64>::zeros(rows, n_params);
    match order {
        0 => {
            for row in 0..n_gamma {
                derivative[(row, row + n_unpenalized)] = 1.0;
            }
        }
        1 => {
            for row in 0..(n_gamma - 1) {
                let delta = validated_log_step(tau, row)?;
                derivative[(row, row + n_unpenalized)] = -1.0 / delta;
                derivative[(row, row + n_unpenalized + 1)] = 1.0 / delta;
            }
        }
        2 => {
            for row in 0..(n_gamma - 2) {
                let delta = validated_log_step(tau, row)?;
                let scale = if row == 0 || row + 1 == n_gamma - 2 {
                    2.0 / (delta * delta)
                } else {
                    1.0 / (delta * delta)
                };
                derivative[(row, row + n_unpenalized)] = scale;
                derivative[(row, row + n_unpenalized + 1)] = -2.0 * scale;
                derivative[(row, row + n_unpenalized + 2)] = scale;
            }
        }
        _ => unreachable!(),
    }
    Ok(derivative.transpose() * derivative)
}

fn validated_log_step(tau: &[f64], row: usize) -> Result<f64> {
    let delta = (tau[row + 1] / tau[row]).ln();
    if delta <= 0.0 || !delta.is_finite() {
        bail!("tau grid must be strictly increasing for DRT regularization");
    }
    Ok(delta)
}
