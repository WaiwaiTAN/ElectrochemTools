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
    if let Some(solution) = system.clone().lu().solve(&atb) {
        return Ok(solution);
    }

    let jitter = (system.trace().abs() / system.nrows().max(1) as f64).max(1.0) * 1.0e-12;
    let stabilized = system + DMatrix::<f64>::identity(a.ncols(), a.ncols()).scale(jitter);
    stabilized
        .lu()
        .solve(&atb)
        .context("Tikhonov system is singular even after numerical stabilization")
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
    let penalty = penalty_matrix(n_gamma, order)?;
    let h = a.transpose() * a + penalty.scale(lambda);
    let g = a.transpose() * b;
    let mut x = solve_tikhonov(a, b, lambda, n_gamma, order)?;
    for idx in 1..x.len() {
        x[idx] = x[idx].max(0.0);
    }

    let step = 1.0 / estimate_largest_eigenvalue(&h).max(1.0e-12);
    let mut previous_objective = objective(&h, &g, &x);
    for _ in 0..max_iter {
        let gradient = &h * &x - &g;
        let mut next = &x - gradient.scale(step);
        for idx in 1..next.len() {
            next[idx] = next[idx].max(0.0);
        }
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
    if a.ncols() != n_gamma + 1 {
        bail!("A column count must be 1 + number of gamma values");
    }
    let penalty = penalty_matrix(n_gamma, order)?;
    let h = a.transpose() * a + penalty.scale(lambda);
    let g = a.transpose() * b;
    let n_params = n_gamma + 1;
    let mut x = DVector::<f64>::zeros(n_params);
    let mut free = vec![false; n_params];
    free[0] = true;

    let unconstrained = solve_tikhonov(a, b, lambda, n_gamma, order)?;
    x[0] = unconstrained[0];
    for idx in 1..n_params {
        if unconstrained[idx] > tol {
            free[idx] = true;
            x[idx] = unconstrained[idx];
        }
    }

    for _ in 0..max_iter {
        let free_indices: Vec<usize> = free
            .iter()
            .enumerate()
            .filter_map(|(idx, is_free)| (*is_free).then_some(idx))
            .collect();
        let candidate = solve_free_subproblem(&h, &g, &free_indices, n_params)?;

        if let Some(blocking_idx) = most_negative_free_gamma(&candidate, &free, tol) {
            let alpha = step_to_boundary(&x, &candidate, &free).unwrap_or(1.0);
            x = &x + (&candidate - &x).scale(alpha);
            for idx in 1..n_params {
                if x[idx] <= tol {
                    x[idx] = 0.0;
                    free[idx] = false;
                }
            }
            free[blocking_idx] = false;
            continue;
        }

        x = candidate;
        for idx in 1..n_params {
            if x[idx] <= tol {
                x[idx] = 0.0;
                free[idx] = false;
            }
        }

        let gradient = &h * &x - &g;
        let mut best_idx = None;
        let mut most_negative_gradient = -tol;
        for idx in 1..n_params {
            if !free[idx] && gradient[idx] < most_negative_gradient {
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
    let solution_free = solve_spd_or_lu(&h_free, &g_free)?;
    let mut x = DVector::<f64>::zeros(n_params);
    for (local, &global) in free_indices.iter().enumerate() {
        x[global] = solution_free[local];
    }
    Ok(x)
}

fn solve_spd_or_lu(system: &DMatrix<f64>, rhs: &DVector<f64>) -> Result<DVector<f64>> {
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
        .context("active-set subproblem is singular even after numerical stabilization")
}

fn most_negative_free_gamma(x: &DVector<f64>, free: &[bool], tol: f64) -> Option<usize> {
    (1..x.len())
        .filter(|&idx| free[idx] && x[idx] < -tol)
        .min_by(|&a, &b| x[a].total_cmp(&x[b]))
}

fn step_to_boundary(
    current: &DVector<f64>,
    candidate: &DVector<f64>,
    free: &[bool],
) -> Option<f64> {
    (1..current.len())
        .filter(|&idx| free[idx] && candidate[idx] <= 0.0)
        .filter_map(|idx| {
            let denom = current[idx] - candidate[idx];
            (denom > 0.0).then_some(current[idx] / denom)
        })
        .min_by(|a, b| a.total_cmp(b))
        .map(|alpha| alpha.clamp(0.0, 1.0))
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
