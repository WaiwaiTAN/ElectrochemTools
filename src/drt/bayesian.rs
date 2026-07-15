use anyhow::{Result, bail};
use nalgebra::{DMatrix, DVector};
use rand::{SeedableRng, rngs::StdRng};
use rand_distr::{Distribution, StandardNormal};
use serde::Serialize;
use std::f64::consts::{FRAC_PI_2, TAU};

const LOWER_PROBABILITY: f64 = 0.005;
const UPPER_PROBABILITY: f64 = 0.995;
const ROOT_EPSILON: f64 = 1.0e-12;
const REHIT_EPSILON: f64 = 1.0e-9;
const FEASIBILITY_TOLERANCE: f64 = 1.0e-10;

#[derive(Debug, Clone, Copy, Serialize)]
pub struct BayesianSettings {
    /// Number of HMC states per chain, including the initial state and burn-in.
    pub samples: usize,
    pub burn_in: usize,
    pub seed: u64,
    pub chains: usize,
}

impl Default for BayesianSettings {
    fn default() -> Self {
        Self {
            samples: 5_000,
            burn_in: 500,
            seed: 0,
            chains: 1,
        }
    }
}

impl BayesianSettings {
    pub fn validate(self) -> Result<()> {
        validate_settings(self)
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct BayesianSamplingResult {
    pub chains: usize,
    pub samples_per_chain: usize,
    pub burn_in: usize,
    pub retained_samples_per_chain: usize,
    pub total_samples: usize,
    pub total_retained_samples: usize,
    pub seed: u64,
    pub chain_seeds: Vec<u64>,
    pub lower_probability: f64,
    pub upper_probability: f64,
    pub noise_std: f64,
    pub bounce_count: usize,
    pub coefficient_r_hat: Vec<f64>,
    pub coefficient_effective_sample_size: Vec<f64>,
    pub max_r_hat: f64,
    pub min_effective_sample_size: f64,
    pub diagnostics_qualified: bool,
    pub coefficient_mean: Vec<f64>,
    pub coefficient_lower: Vec<f64>,
    pub coefficient_upper: Vec<f64>,
}

pub(super) fn sample_nonnegative_posterior(
    matrix: &DMatrix<f64>,
    observations: &DVector<f64>,
    ridge_coefficients: &DVector<f64>,
    lambda: f64,
    penalty: &DMatrix<f64>,
    n_unpenalized: usize,
    settings: BayesianSettings,
) -> Result<BayesianSamplingResult> {
    validate_settings(settings)?;
    if matrix.nrows() != observations.len()
        || matrix.ncols() != ridge_coefficients.len()
        || penalty.shape() != (matrix.ncols(), matrix.ncols())
        || n_unpenalized >= matrix.ncols()
    {
        bail!("invalid matrix dimensions for Bayesian DRT sampling");
    }

    let residual = matrix * ridge_coefficients - observations;
    let residual_mean = residual.iter().sum::<f64>() / residual.len() as f64;
    let residual_ss = residual
        .iter()
        .map(|value| (value - residual_mean).powi(2))
        .sum::<f64>();
    let noise_variance = residual_ss / (residual.len().saturating_sub(1).max(1) as f64);
    if !noise_variance.is_finite() || noise_variance <= f64::EPSILON {
        bail!("Bayesian DRT requires a finite, non-zero residual noise estimate");
    }

    // This is the posterior construction used by DRTtools: the observation
    // variance scales both the likelihood and Tikhonov prior, and only the
    // gamma block is passed to the nonnegative sampler.
    let precision_unscaled = matrix.transpose() * matrix + penalty.scale(lambda);
    let rhs_unscaled = matrix.transpose() * observations;
    let posterior_mean_full = precision_unscaled
        .clone()
        .cholesky()
        .ok_or_else(|| {
            anyhow::anyhow!("Bayesian DRT posterior precision is not positive definite")
        })?
        .solve(&rhs_unscaled);

    let n_gamma = matrix.ncols() - n_unpenalized;
    let gamma_precision = precision_unscaled
        .view((n_unpenalized, n_unpenalized), (n_gamma, n_gamma))
        .into_owned()
        .scale(1.0 / noise_variance);
    let gamma_covariance = gamma_precision
        .try_inverse()
        .ok_or_else(|| anyhow::anyhow!("Bayesian DRT gamma precision is singular"))?;
    let gamma_covariance = 0.5 * (&gamma_covariance + gamma_covariance.transpose());
    let covariance_factor = gamma_covariance
        .cholesky()
        .ok_or_else(|| anyhow::anyhow!("Bayesian DRT gamma covariance is not positive definite"))?
        .l();
    let gamma_mean = posterior_mean_full
        .rows(n_unpenalized, n_gamma)
        .into_owned();

    let scale = ridge_coefficients
        .rows(n_unpenalized, n_gamma)
        .iter()
        .copied()
        .fold(1.0_f64, |current, value| current.max(value.abs()));
    let positive_floor = 100.0 * f64::EPSILON * scale;
    let initial = DVector::from_iterator(
        n_gamma,
        ridge_coefficients
            .iter()
            .skip(n_unpenalized)
            .map(|value| value.max(positive_floor)),
    );
    let initial_white = covariance_factor
        .clone()
        .lu()
        .solve(&(initial - &gamma_mean))
        .ok_or_else(|| anyhow::anyhow!("failed to whiten the Bayesian DRT initial state"))?;
    let transform_matrix = covariance_factor;
    let transform_offset = gamma_mean;
    let mut constraint_matrix = transform_matrix.clone();
    let mut constraint_offset = transform_offset.clone();
    // Scaling a half-space inequality by a positive constant does not alter
    // the feasible region or reflection. Unit row norms avoid catastrophic
    // cancellation for ill-conditioned Gaussian RBF covariance factors.
    for row in 0..constraint_matrix.nrows() {
        let norm = constraint_matrix.row(row).norm();
        if !norm.is_finite() || norm <= f64::MIN_POSITIVE {
            bail!("Bayesian DRT sampler encountered a degenerate constraint");
        }
        constraint_matrix.row_mut(row).scale_mut(1.0 / norm);
        constraint_offset[row] /= norm;
    }

    let chain_seeds = (0..settings.chains)
        .map(|chain| derive_chain_seed(settings.seed, chain))
        .collect::<Vec<_>>();
    let chain_outputs = if settings.chains == 1 {
        let output = exact_hmc(
            &constraint_matrix,
            &constraint_offset,
            initial_white,
            settings.samples,
            chain_seeds[0],
        )?;
        vec![output]
    } else {
        std::thread::scope(|scope| -> Result<Vec<_>> {
            let mut handles = Vec::with_capacity(settings.chains);
            for &chain_seed in &chain_seeds {
                let constraints = &constraint_matrix;
                let offset = &constraint_offset;
                let initial = initial_white.clone();
                handles.push(scope.spawn(move || {
                    exact_hmc(constraints, offset, initial, settings.samples, chain_seed)
                }));
            }
            handles
                .into_iter()
                .map(|handle| {
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("Bayesian DRT worker chain panicked"))?
                })
                .collect()
        })?
    };
    let bounce_count = chain_outputs.iter().map(|(_, count)| count).sum();
    let retained_chains = chain_outputs
        .iter()
        .map(|(states, _)| {
            states[settings.burn_in..]
                .iter()
                .map(|state| {
                    DVector::from_fn(n_gamma, |coefficient, _| {
                        (transform_matrix.row(coefficient) * state)[0]
                            + transform_offset[coefficient]
                    })
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut coefficient_mean = vec![0.0; n_gamma];
    let mut coefficient_lower = vec![0.0; n_gamma];
    let mut coefficient_upper = vec![0.0; n_gamma];
    let mut coefficient_r_hat = vec![0.0; n_gamma];
    let mut coefficient_effective_sample_size = vec![0.0; n_gamma];
    for coefficient in 0..n_gamma {
        let coefficient_chains = retained_chains
            .iter()
            .map(|chain| {
                chain
                    .iter()
                    .map(|state| state[coefficient])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut values = coefficient_chains
            .iter()
            .flatten()
            .copied()
            .collect::<Vec<_>>();
        coefficient_mean[coefficient] = values.iter().sum::<f64>() / values.len() as f64;
        values.sort_by(f64::total_cmp);
        coefficient_lower[coefficient] = quantile_r5(&values, LOWER_PROBABILITY);
        coefficient_upper[coefficient] = quantile_r5(&values, UPPER_PROBABILITY);
        coefficient_r_hat[coefficient] = split_r_hat(&coefficient_chains);
        coefficient_effective_sample_size[coefficient] =
            effective_sample_size(&coefficient_chains, coefficient_r_hat[coefficient]);
    }
    let max_r_hat = coefficient_r_hat
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let min_effective_sample_size = coefficient_effective_sample_size
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let diagnostics_qualified = max_r_hat <= 1.01 && min_effective_sample_size >= 400.0;

    Ok(BayesianSamplingResult {
        chains: settings.chains,
        samples_per_chain: settings.samples,
        burn_in: settings.burn_in,
        retained_samples_per_chain: settings.samples - settings.burn_in,
        total_samples: settings.chains * settings.samples,
        total_retained_samples: settings.chains * (settings.samples - settings.burn_in),
        seed: settings.seed,
        chain_seeds,
        lower_probability: LOWER_PROBABILITY,
        upper_probability: UPPER_PROBABILITY,
        noise_std: noise_variance.sqrt(),
        bounce_count,
        coefficient_r_hat,
        coefficient_effective_sample_size,
        max_r_hat,
        min_effective_sample_size,
        diagnostics_qualified,
        coefficient_mean,
        coefficient_lower,
        coefficient_upper,
    })
}

fn validate_settings(settings: BayesianSettings) -> Result<()> {
    if settings.samples < 1_000 {
        bail!("Bayesian DRT requires at least 1000 samples per chain");
    }
    if settings.burn_in >= settings.samples {
        bail!("Bayesian DRT burn-in must be smaller than the total sample count");
    }
    if settings.samples - settings.burn_in < 2 {
        bail!("Bayesian DRT requires at least two retained samples");
    }
    if settings.chains == 0 || settings.chains > 64 {
        bail!("Bayesian DRT chain count must be between 1 and 64");
    }
    settings
        .chains
        .checked_mul(settings.samples)
        .ok_or_else(|| anyhow::anyhow!("Bayesian DRT total sample count overflows usize"))?;
    Ok(())
}

fn derive_chain_seed(base_seed: u64, chain: usize) -> u64 {
    if chain == 0 {
        return base_seed;
    }
    let mut value = base_seed.wrapping_add((chain as u64).wrapping_mul(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

fn split_r_hat(chains: &[Vec<f64>]) -> f64 {
    let half_length = chains.iter().map(Vec::len).min().unwrap_or(0) / 2;
    if half_length < 2 {
        return f64::INFINITY;
    }
    let mut means = Vec::with_capacity(2 * chains.len());
    let mut variances = Vec::with_capacity(2 * chains.len());
    for chain in chains {
        for split in [&chain[..half_length], &chain[chain.len() - half_length..]] {
            let mean = split.iter().sum::<f64>() / half_length as f64;
            let variance = split
                .iter()
                .map(|value| (value - mean).powi(2))
                .sum::<f64>()
                / (half_length - 1) as f64;
            means.push(mean);
            variances.push(variance);
        }
    }
    let within = variances.iter().sum::<f64>() / variances.len() as f64;
    let mean_of_means = means.iter().sum::<f64>() / means.len() as f64;
    let between = half_length as f64
        * means
            .iter()
            .map(|mean| (mean - mean_of_means).powi(2))
            .sum::<f64>()
        / (means.len() - 1) as f64;
    if within <= f64::EPSILON {
        return if between <= f64::EPSILON {
            1.0
        } else {
            f64::INFINITY
        };
    }
    let variance_estimate = ((half_length - 1) as f64 * within + between) / half_length as f64;
    (variance_estimate / within).max(1.0).sqrt()
}

fn effective_sample_size(chains: &[Vec<f64>], r_hat: f64) -> f64 {
    let chain_length = chains.iter().map(Vec::len).min().unwrap_or(0);
    if chain_length < 3 || chains.is_empty() {
        return 0.0;
    }
    let means = chains
        .iter()
        .map(|chain| chain[..chain_length].iter().sum::<f64>() / chain_length as f64)
        .collect::<Vec<_>>();
    let variance = chains
        .iter()
        .zip(&means)
        .map(|(chain, mean)| {
            chain[..chain_length]
                .iter()
                .map(|value| (value - mean).powi(2))
                .sum::<f64>()
                / (chain_length - 1) as f64
        })
        .sum::<f64>()
        / chains.len() as f64;
    let total = (chains.len() * chain_length) as f64;
    if variance <= f64::EPSILON {
        return total;
    }

    let autocorrelation = |lag: usize| {
        chains
            .iter()
            .zip(&means)
            .map(|(chain, mean)| {
                (0..(chain_length - lag))
                    .map(|index| (chain[index] - mean) * (chain[index + lag] - mean))
                    .sum::<f64>()
                    / (chain_length - lag) as f64
            })
            .sum::<f64>()
            / chains.len() as f64
            / variance
    };

    let mut paired_sum = 0.0;
    let max_lag = (chain_length - 1).min(1_000);
    let mut lag = 1;
    while lag <= max_lag {
        let pair = autocorrelation(lag)
            + if lag < max_lag {
                autocorrelation(lag + 1)
            } else {
                0.0
            };
        if pair <= 0.0 {
            break;
        }
        paired_sum += pair;
        lag += 2;
    }
    let r_hat_penalty = if r_hat.is_finite() {
        r_hat * r_hat
    } else {
        f64::INFINITY
    };
    (total / ((1.0 + 2.0 * paired_sum) * r_hat_penalty)).clamp(0.0, total)
}

/// Exact HMC for a standard normal subject to `constraints * x + offset >= 0`.
/// The harmonic trajectory and specular wall reflections follow Pakman and
/// Paninski (2014), matching the sampler bundled with DRTtools.
fn exact_hmc(
    constraints: &DMatrix<f64>,
    offset: &DVector<f64>,
    initial: DVector<f64>,
    samples: usize,
    seed: u64,
) -> Result<(Vec<DVector<f64>>, usize)> {
    if constraints.nrows() != offset.len() || constraints.ncols() != initial.len() {
        bail!("invalid exact-HMC constraint dimensions");
    }
    if (constraints * &initial + offset)
        .iter()
        .any(|value| *value < -FEASIBILITY_TOLERANCE)
    {
        bail!("Bayesian DRT sampler initial state is infeasible");
    }

    let row_norm_squared = (0..constraints.nrows())
        .map(|row| constraints.row(row).norm_squared())
        .collect::<Vec<_>>();
    if row_norm_squared.iter().any(|value| *value <= f64::EPSILON) {
        bail!("Bayesian DRT sampler encountered a degenerate constraint");
    }

    let mut rng = StdRng::seed_from_u64(seed);
    let normal = StandardNormal;
    let mut current = initial;
    let mut states = Vec::with_capacity(samples);
    states.push(current.clone());
    let mut bounce_count = 0usize;
    let mut rejected_trajectories = 0usize;
    let mut largest_feasibility_violation = 0.0_f64;

    while states.len() < samples {
        let mut velocity = DVector::from_fn(current.len(), |_, _| normal.sample(&mut rng));
        let mut position = current.clone();
        let mut remaining = FRAC_PI_2;
        let mut reflections = 0usize;
        let mut last_hit = None;

        while remaining > ROOT_EPSILON {
            let mut hit_time = remaining + 1.0;
            let mut hit_row = None;
            for row in 0..constraints.nrows() {
                let normal_row = constraints.row(row);
                let sin_coefficient = (normal_row * &velocity)[0];
                let cos_coefficient = (normal_row * &position)[0];
                if let Some(time) = first_wall_time(sin_coefficient, cos_coefficient, offset[row])
                    && !(last_hit == Some(row) && time <= REHIT_EPSILON)
                    && time < hit_time
                    && time <= remaining + ROOT_EPSILON
                {
                    hit_time = time;
                    hit_row = Some(row);
                }
            }

            let movement = hit_row.map_or(remaining, |_| hit_time.min(remaining));
            let old_position = position;
            position = &velocity * movement.sin() + &old_position * movement.cos();
            velocity = velocity * movement.cos() - old_position * movement.sin();
            remaining = (remaining - movement).max(0.0);

            let segment_feasibility = constraints * &position + offset;
            if let Some((bad_row, bad_value)) = segment_feasibility
                .iter()
                .enumerate()
                .min_by(|(_, left), (_, right)| left.total_cmp(right))
                .filter(|(_, value)| **value < -1.0e-6)
            {
                bail!(
                    "Bayesian DRT HMC crossed constraint {bad_row} by {:.3e} (selected wall {:?}, movement {:.3e})",
                    -*bad_value,
                    hit_row,
                    movement
                );
            }

            let Some(row) = hit_row.filter(|_| hit_time <= movement + ROOT_EPSILON) else {
                break;
            };
            let normal_row = constraints.row(row);
            // Ill-conditioned Gaussian RBF covariance factors can accumulate
            // cancellation error at a wall. Put the position back on the
            // analytically hit plane before reflecting the velocity.
            let wall_error = (normal_row * &position)[0] + offset[row];
            position -= normal_row.transpose() * (wall_error / row_norm_squared[row]);
            let projection = (normal_row * &velocity)[0] / row_norm_squared[row];
            velocity -= normal_row.transpose() * (2.0 * projection);
            last_hit = Some(row);
            bounce_count += 1;
            reflections += 1;
            if reflections > 100_000 {
                bail!("Bayesian DRT HMC trajectory exceeded the reflection limit");
            }
        }

        let feasibility = constraints * &position + offset;
        if feasibility
            .iter()
            .all(|value| *value >= -FEASIBILITY_TOLERANCE)
        {
            current = position;
            states.push(current.clone());
        } else {
            largest_feasibility_violation = largest_feasibility_violation.max(
                feasibility
                    .iter()
                    .map(|value| (-value).max(0.0))
                    .fold(0.0, f64::max),
            );
            rejected_trajectories += 1;
            if rejected_trajectories > 10_000 {
                bail!(
                    "Bayesian DRT HMC exceeded the rejected-trajectory limit (largest feasibility violation {:.3e})",
                    largest_feasibility_violation
                );
            }
        }
    }
    Ok((states, bounce_count))
}

fn first_wall_time(sin_coefficient: f64, cos_coefficient: f64, offset: f64) -> Option<f64> {
    if cos_coefficient + offset <= FEASIBILITY_TOLERANCE && sin_coefficient < 0.0 {
        return Some(0.0);
    }
    let amplitude = sin_coefficient.hypot(cos_coefficient);
    if amplitude <= f64::EPSILON || offset.abs() > amplitude {
        return None;
    }
    let phase = (-sin_coefficient).atan2(cos_coefficient);
    let angle = (-offset / amplitude).clamp(-1.0, 1.0).acos();
    let root = (-phase + angle).rem_euclid(TAU);
    Some(if root <= ROOT_EPSILON {
        root + TAU
    } else {
        root
    })
}

fn quantile_r5(sorted: &[f64], probability: f64) -> f64 {
    debug_assert!(sorted.len() >= 2);
    let n = sorted.len() as f64;
    if probability < 0.5 / n {
        return sorted[0];
    }
    if probability >= (n - 0.5) / n {
        return sorted[sorted.len() - 1];
    }
    let h = n * probability + 0.5;
    let lower_one_based = h.floor() as usize;
    let fraction = h - h.floor();
    let lower = sorted[lower_one_based - 1];
    let upper = sorted[lower_one_based];
    lower + fraction * (upper - lower)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r5_quantile_interpolates_using_matlab_convention() {
        let values = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(quantile_r5(&values, 0.5), 2.5);
        assert_eq!(quantile_r5(&values, 0.0), 1.0);
        assert_eq!(quantile_r5(&values, 1.0), 4.0);
    }

    #[test]
    fn exact_hmc_is_seeded_and_respects_positive_constraints() {
        let constraints = DMatrix::identity(2, 2);
        let offset = DVector::zeros(2);
        let initial = DVector::from_vec(vec![0.5, 0.5]);
        let (first, _) = exact_hmc(&constraints, &offset, initial.clone(), 20, 42).unwrap();
        let (second, _) = exact_hmc(&constraints, &offset, initial, 20, 42).unwrap();
        assert_eq!(first, second);
        assert!(first.iter().flatten().all(|value| *value >= -1.0e-10));
    }

    #[test]
    fn exact_hmc_respects_correlated_offset_constraints() {
        let constraints =
            DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 0.0, 0.95, 0.1, 0.0, -0.7, 0.4, 0.05]);
        let offset = DVector::from_vec(vec![-1.0, -0.8, -0.6]);
        let feasible = DVector::from_vec(vec![0.2, 0.3, 0.4]);
        let initial = constraints
            .clone()
            .lu()
            .solve(&(feasible - &offset))
            .unwrap();
        let (states, _) = exact_hmc(&constraints, &offset, initial, 1_000, 42).unwrap();
        assert!(states.iter().all(|state| {
            (&constraints * state + &offset)
                .iter()
                .all(|value| *value >= -1.0e-10)
        }));
    }

    #[test]
    fn one_dimensional_sampler_matches_half_normal_mean() {
        let constraints = DMatrix::identity(1, 1);
        let offset = DVector::zeros(1);
        let initial = DVector::from_element(1, 0.5);
        let (states, _) = exact_hmc(&constraints, &offset, initial, 10_000, 7).unwrap();
        let mean = states.iter().skip(100).map(|state| state[0]).sum::<f64>()
            / (states.len() - 100) as f64;
        let expected = (2.0 / std::f64::consts::PI).sqrt();
        assert!((mean - expected).abs() < 0.03, "{mean} vs {expected}");
    }

    #[test]
    fn settings_reject_short_or_empty_retained_chains() {
        assert!(
            BayesianSettings {
                samples: 999,
                burn_in: 100,
                seed: 0,
                chains: 1,
            }
            .validate()
            .is_err()
        );
        assert!(
            BayesianSettings {
                samples: 1_000,
                burn_in: 1_000,
                seed: 0,
                chains: 1,
            }
            .validate()
            .is_err()
        );
    }

    #[test]
    fn chain_seeds_are_stable_distinct_and_keep_first_chain_compatible() {
        let seeds = (0..4)
            .map(|chain| derive_chain_seed(42, chain))
            .collect::<Vec<_>>();
        assert_eq!(seeds[0], 42);
        assert_eq!(
            seeds,
            vec![
                42,
                13679457532755275413,
                2949826092126892291,
                5139283748462763858,
            ]
        );
    }

    #[test]
    fn diagnostics_distinguish_mixed_and_separated_chains() {
        let mixed = (0..4)
            .map(|chain| {
                (0..1_000)
                    .map(|index| ((index * 37 + chain * 101) as f64).sin())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mixed_r_hat = split_r_hat(&mixed);
        assert!(mixed_r_hat < 1.01, "{mixed_r_hat}");
        assert!(effective_sample_size(&mixed, mixed_r_hat) > 400.0);

        let separated = (0..4)
            .map(|chain| {
                (0..1_000)
                    .map(|index| chain as f64 * 10.0 + (index as f64).sin())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        assert!(split_r_hat(&separated) > 1.1);
    }
}
