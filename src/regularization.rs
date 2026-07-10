use anyhow::{Context, Result, bail};
use nalgebra::{DMatrix, DVector};

pub fn solve_tikhonov(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    n_gamma: usize,
    order: usize,
) -> Result<DVector<f64>> {
    solve_tikhonov_general(a, b, lambda, n_gamma, order, 1)
}

pub fn solve_tikhonov_general(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    n_gamma: usize,
    order: usize,
    n_unpenalized: usize,
) -> Result<DVector<f64>> {
    validate_inputs(a, b, lambda, n_gamma, n_unpenalized)?;
    let system = a.transpose() * a
        + penalty_matrix_with_unpenalized(n_gamma, order, n_unpenalized)?.scale(lambda);
    let rhs = a.transpose() * b;
    solve_spd_or_lu(&system, &rhs, "Tikhonov system is singular")
}

pub fn solve_tikhonov_general_with_penalty(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    penalty: &DMatrix<f64>,
) -> Result<DVector<f64>> {
    validate_penalty_inputs(a, b, lambda, penalty)?;
    let system = a.transpose() * a + penalty.scale(lambda);
    let rhs = a.transpose() * b;
    solve_spd_or_lu(&system, &rhs, "Tikhonov system is singular")
}

pub fn solve_tikhonov_projected_nonnegative(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    n_gamma: usize,
    order: usize,
    max_iter: usize,
    tol: f64,
) -> Result<DVector<f64>> {
    let lower_bound_zero = (0..(n_gamma + 1)).map(|idx| idx >= 1).collect::<Vec<_>>();
    solve_tikhonov_projected_with_bounds(
        a,
        b,
        lambda,
        n_gamma,
        order,
        1,
        &lower_bound_zero,
        max_iter,
        tol,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn solve_tikhonov_projected_with_bounds(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    n_gamma: usize,
    order: usize,
    n_unpenalized: usize,
    lower_bound_zero: &[bool],
    max_iter: usize,
    tol: f64,
) -> Result<DVector<f64>> {
    validate_inputs(a, b, lambda, n_gamma, n_unpenalized)?;
    validate_lower_bounds(a.ncols(), lower_bound_zero)?;
    let h = a.transpose() * a
        + penalty_matrix_with_unpenalized(n_gamma, order, n_unpenalized)?.scale(lambda);
    let g = a.transpose() * b;
    let mut x = solve_tikhonov_general(a, b, lambda, n_gamma, order, n_unpenalized)?;
    project_bounds(&mut x, lower_bound_zero);

    let step = 1.0 / estimate_largest_eigenvalue(&h).max(1.0e-12);
    let mut previous_objective = objective(&h, &g, &x);
    for _ in 0..max_iter {
        let gradient = &h * &x - &g;
        let mut next = &x - gradient.scale(step);
        project_bounds(&mut next, lower_bound_zero);
        let next_objective = objective(&h, &g, &next);
        let improvement = (previous_objective - next_objective).abs();
        x = next;
        if improvement <= tol * previous_objective.abs().max(1.0) {
            break;
        }
        previous_objective = next_objective;
    }
    Ok(x)
}

pub fn solve_tikhonov_active_set_nonnegative(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    n_gamma: usize,
    order: usize,
    max_iter: usize,
    tol: f64,
) -> Result<DVector<f64>> {
    let lower_bound_zero = (0..(n_gamma + 1)).map(|idx| idx >= 1).collect::<Vec<_>>();
    solve_tikhonov_active_set_with_bounds(
        a,
        b,
        lambda,
        n_gamma,
        order,
        1,
        &lower_bound_zero,
        max_iter,
        tol,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn solve_tikhonov_active_set_with_bounds(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    n_gamma: usize,
    order: usize,
    n_unpenalized: usize,
    lower_bound_zero: &[bool],
    max_iter: usize,
    tol: f64,
) -> Result<DVector<f64>> {
    validate_inputs(a, b, lambda, n_gamma, n_unpenalized)?;
    let penalty = penalty_matrix_with_unpenalized(n_gamma, order, n_unpenalized)?;
    solve_tikhonov_active_set_with_penalty(a, b, lambda, &penalty, lower_bound_zero, max_iter, tol)
}

pub fn solve_tikhonov_active_set_with_penalty(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    penalty: &DMatrix<f64>,
    lower_bound_zero: &[bool],
    max_iter: usize,
    tol: f64,
) -> Result<DVector<f64>> {
    validate_penalty_inputs(a, b, lambda, penalty)?;
    validate_lower_bounds(a.ncols(), lower_bound_zero)?;
    let h = a.transpose() * a + penalty.scale(lambda);
    let g = a.transpose() * b;
    let n_params = a.ncols();
    let mut x = DVector::<f64>::zeros(n_params);
    let mut free = vec![false; n_params];

    let unconstrained = solve_tikhonov_general_with_penalty(a, b, lambda, penalty)?;
    for idx in 0..n_params {
        if lower_bound_zero[idx] {
            if unconstrained[idx] > tol {
                free[idx] = true;
                x[idx] = unconstrained[idx];
            }
        } else {
            free[idx] = true;
            x[idx] = unconstrained[idx];
        }
    }

    for _ in 0..max_iter {
        let free_indices = free
            .iter()
            .enumerate()
            .filter_map(|(idx, is_free)| (*is_free).then_some(idx))
            .collect::<Vec<_>>();
        let candidate = solve_free_subproblem(&h, &g, &free_indices, n_params)?;

        if let Some(blocking_idx) =
            most_negative_free_bounded(&candidate, &free, lower_bound_zero, tol)
        {
            let alpha = step_to_boundary(&x, &candidate, &free, lower_bound_zero).unwrap_or(1.0);
            x = &x + (&candidate - &x).scale(alpha);
            project_inactive(&mut x, &mut free, lower_bound_zero, tol);
            free[blocking_idx] = false;
            continue;
        }

        x = candidate;
        project_inactive(&mut x, &mut free, lower_bound_zero, tol);

        let gradient = &h * &x - &g;
        let mut best_idx = None;
        let mut most_negative_gradient = -tol;
        for idx in 0..n_params {
            if lower_bound_zero[idx] && !free[idx] && gradient[idx] < most_negative_gradient {
                most_negative_gradient = gradient[idx];
                best_idx = Some(idx);
            }
        }

        if let Some(idx) = best_idx {
            free[idx] = true;
        } else {
            return Ok(x);
        }
    }

    Ok(x)
}

pub fn penalty_matrix(n_gamma: usize, order: usize) -> Result<DMatrix<f64>> {
    penalty_matrix_with_unpenalized(n_gamma, order, 1)
}

pub fn penalty_matrix_with_unpenalized(
    n_gamma: usize,
    order: usize,
    n_unpenalized: usize,
) -> Result<DMatrix<f64>> {
    let n_params = n_gamma + n_unpenalized;
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
                l[(row, row + n_unpenalized)] = 1.0;
            }
        }
        1 => {
            for row in 0..(n_gamma - 1) {
                l[(row, row + n_unpenalized)] = -1.0;
                l[(row, row + n_unpenalized + 1)] = 1.0;
            }
        }
        2 => {
            for row in 0..(n_gamma - 2) {
                l[(row, row + n_unpenalized)] = 1.0;
                l[(row, row + n_unpenalized + 1)] = -2.0;
                l[(row, row + n_unpenalized + 2)] = 1.0;
            }
        }
        _ => unreachable!(),
    }
    Ok(l.transpose() * l)
}

fn validate_inputs(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    n_gamma: usize,
    n_unpenalized: usize,
) -> Result<()> {
    if a.nrows() != b.len() {
        bail!("A row count must match b length");
    }
    if a.ncols() != n_gamma + n_unpenalized {
        bail!("A column count must be number of gamma values plus unpenalized variables");
    }
    if lambda < 0.0 || !lambda.is_finite() {
        bail!("lambda must be finite and non-negative");
    }
    Ok(())
}

fn validate_lower_bounds(n_params: usize, lower_bound_zero: &[bool]) -> Result<()> {
    if lower_bound_zero.len() != n_params {
        bail!("lower bound vector length must match parameter count");
    }
    Ok(())
}

fn validate_penalty_inputs(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    penalty: &DMatrix<f64>,
) -> Result<()> {
    if a.nrows() != b.len() {
        bail!("A row count must match b length");
    }
    if penalty.nrows() != a.ncols() || penalty.ncols() != a.ncols() {
        bail!("penalty matrix must be square with one row/column per parameter");
    }
    if lambda < 0.0 || !lambda.is_finite() {
        bail!("lambda must be finite and non-negative");
    }
    Ok(())
}

fn solve_free_subproblem(
    h: &DMatrix<f64>,
    g: &DVector<f64>,
    free_indices: &[usize],
    n_params: usize,
) -> Result<DVector<f64>> {
    let mut h_free = DMatrix::<f64>::zeros(free_indices.len(), free_indices.len());
    let mut g_free = DVector::<f64>::zeros(free_indices.len());
    for (row_local, &row_global) in free_indices.iter().enumerate() {
        g_free[row_local] = g[row_global];
        for (col_local, &col_global) in free_indices.iter().enumerate() {
            h_free[(row_local, col_local)] = h[(row_global, col_global)];
        }
    }
    let solution_free = solve_spd_or_lu(
        &h_free,
        &g_free,
        "active-set subproblem is singular even after numerical stabilization",
    )?;
    let mut x = DVector::<f64>::zeros(n_params);
    for (local, &global) in free_indices.iter().enumerate() {
        x[global] = solution_free[local];
    }
    Ok(x)
}

fn solve_spd_or_lu(
    system: &DMatrix<f64>,
    rhs: &DVector<f64>,
    singular_message: &str,
) -> Result<DVector<f64>> {
    if let Some(cholesky) = system.clone().cholesky() {
        return Ok(cholesky.solve(rhs));
    }
    if let Some(solution) = system.clone().lu().solve(rhs) {
        return Ok(solution);
    }
    let jitter = (system.trace().abs() / system.nrows().max(1) as f64).max(1.0) * 1.0e-12;
    let stabilized =
        system + DMatrix::<f64>::identity(system.nrows(), system.ncols()).scale(jitter);
    stabilized
        .lu()
        .solve(rhs)
        .context(singular_message.to_string())
}

fn most_negative_free_bounded(
    x: &DVector<f64>,
    free: &[bool],
    lower_bound_zero: &[bool],
    tol: f64,
) -> Option<usize> {
    (0..x.len())
        .filter(|&idx| lower_bound_zero[idx] && free[idx] && x[idx] < -tol)
        .min_by(|&a, &b| x[a].total_cmp(&x[b]))
}

fn step_to_boundary(
    current: &DVector<f64>,
    candidate: &DVector<f64>,
    free: &[bool],
    lower_bound_zero: &[bool],
) -> Option<f64> {
    (0..current.len())
        .filter(|&idx| lower_bound_zero[idx] && free[idx] && candidate[idx] <= 0.0)
        .filter_map(|idx| {
            let denom = current[idx] - candidate[idx];
            (denom > 0.0).then_some(current[idx] / denom)
        })
        .min_by(|a, b| a.total_cmp(b))
        .map(|alpha| alpha.clamp(0.0, 1.0))
}

fn project_inactive(x: &mut DVector<f64>, free: &mut [bool], lower_bound_zero: &[bool], tol: f64) {
    for idx in 0..x.len() {
        if lower_bound_zero[idx] && x[idx] <= tol {
            x[idx] = 0.0;
            free[idx] = false;
        }
    }
}

fn project_bounds(x: &mut DVector<f64>, lower_bound_zero: &[bool]) {
    for idx in 0..x.len() {
        if lower_bound_zero[idx] {
            x[idx] = x[idx].max(0.0);
        }
    }
}

fn objective(h: &DMatrix<f64>, g: &DVector<f64>, x: &DVector<f64>) -> f64 {
    0.5 * x.dot(&(h * x)) - g.dot(x)
}

fn estimate_largest_eigenvalue(matrix: &DMatrix<f64>) -> f64 {
    let mut v = DVector::<f64>::from_element(matrix.ncols(), 1.0 / matrix.ncols() as f64);
    let mut eigenvalue = 0.0;
    for _ in 0..64 {
        let next = matrix * &v;
        let norm = next.norm();
        if norm <= 0.0 || !norm.is_finite() {
            return 1.0;
        }
        v = next.scale(1.0 / norm);
        eigenvalue = v.dot(&(matrix * &v));
    }
    eigenvalue.abs()
}
