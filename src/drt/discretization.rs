use super::kernel::{BasisKernelMatrices, piecewise_linear_kernel_matrices};
use super::regularization::piecewise_linear_penalty;
use anyhow::{Result, bail};
use nalgebra::{DMatrix, DVector};
use serde::Serialize;
use std::f64::consts::{LN_2, PI};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, clap::ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum DrtBasis {
    #[default]
    PiecewiseLinear,
    Gaussian,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, clap::ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum ShapeControl {
    #[default]
    Fwhm,
    ShapeFactor,
}

/// Basis-specific DRT discretization operations.
///
/// The solver only consumes the matrices returned by this trait, so adding a
/// basis does not require another constrained or combined Re/Im fitting path.
pub trait DrtDiscretization {
    fn basis(&self) -> DrtBasis;
    fn tau(&self) -> &[f64];
    fn epsilon(&self) -> Option<f64>;
    fn kernel_matrices(&self, frequency_hz: &[f64]) -> Result<BasisKernelMatrices>;
    fn penalty(&self, order: usize, n_unpenalized: usize) -> Result<DMatrix<f64>>;
    fn gamma_mapping_matrix(&self) -> DMatrix<f64>;

    fn map_coefficients_to_gamma(&self, coefficients: &DVector<f64>) -> Result<DVector<f64>> {
        if coefficients.len() != self.tau().len() {
            bail!(
                "basis coefficient count {} does not match tau-center count {}",
                coefficients.len(),
                self.tau().len()
            );
        }
        Ok(self.gamma_mapping_matrix() * coefficients)
    }
}

#[derive(Debug, Clone)]
pub struct PiecewiseLinearDiscretization {
    tau: Vec<f64>,
}

impl PiecewiseLinearDiscretization {
    pub fn new(tau: Vec<f64>) -> Result<Self> {
        validate_tau(&tau)?;
        Ok(Self { tau })
    }
}

impl DrtDiscretization for PiecewiseLinearDiscretization {
    fn basis(&self) -> DrtBasis {
        DrtBasis::PiecewiseLinear
    }

    fn tau(&self) -> &[f64] {
        &self.tau
    }

    fn epsilon(&self) -> Option<f64> {
        None
    }

    fn kernel_matrices(&self, frequency_hz: &[f64]) -> Result<BasisKernelMatrices> {
        Ok(piecewise_linear_kernel_matrices(frequency_hz, &self.tau))
    }

    fn penalty(&self, order: usize, n_unpenalized: usize) -> Result<DMatrix<f64>> {
        piecewise_linear_penalty(&self.tau, order, n_unpenalized)
    }

    fn gamma_mapping_matrix(&self) -> DMatrix<f64> {
        DMatrix::identity(self.tau.len(), self.tau.len())
    }
}

#[derive(Debug, Clone)]
pub struct GaussianDiscretization {
    tau: Vec<f64>,
    epsilon: f64,
}

impl GaussianDiscretization {
    pub fn new(tau: Vec<f64>, shape_control: ShapeControl, shape_coefficient: f64) -> Result<Self> {
        validate_tau(&tau)?;
        if !shape_coefficient.is_finite() || shape_coefficient <= 0.0 {
            bail!("shape coefficient must be finite and strictly positive");
        }
        let epsilon = match shape_control {
            ShapeControl::Fwhm => gaussian_epsilon_fwhm(&tau, shape_coefficient)?,
            ShapeControl::ShapeFactor => shape_coefficient,
        };
        Ok(Self { tau, epsilon })
    }
}

impl DrtDiscretization for GaussianDiscretization {
    fn basis(&self) -> DrtBasis {
        DrtBasis::Gaussian
    }

    fn tau(&self) -> &[f64] {
        &self.tau
    }

    fn epsilon(&self) -> Option<f64> {
        Some(self.epsilon)
    }

    fn kernel_matrices(&self, frequency_hz: &[f64]) -> Result<BasisKernelMatrices> {
        let mut real = DMatrix::zeros(frequency_hz.len(), self.tau.len());
        let mut imaginary = DMatrix::zeros(frequency_hz.len(), self.tau.len());
        for &frequency in frequency_hz {
            if !frequency.is_finite() || frequency <= 0.0 {
                bail!("Gaussian DRT kernels require finite positive frequencies");
            }
        }

        // DRTtools deliberately applies a Toeplitz approximation when the
        // inverse-frequency centers are log-spaced to within 1%. Reproducing
        // that branch matters for data whose printed frequencies are only
        // approximately, rather than exactly, logarithmically spaced.
        if frequency_hz.len() == self.tau.len()
            && uses_drttools_toeplitz_approximation(frequency_hz)
        {
            let mut first_column = Vec::with_capacity(frequency_hz.len());
            let mut first_row = Vec::with_capacity(self.tau.len());
            for &frequency in frequency_hz {
                let (kernel_real, kernel_imaginary_positive) =
                    gaussian_debye_kernel(2.0 * PI * frequency * self.tau[0], self.epsilon);
                first_column.push((kernel_real, -kernel_imaginary_positive));
            }
            for &tau in &self.tau {
                let (kernel_real, kernel_imaginary_positive) =
                    gaussian_debye_kernel(2.0 * PI * frequency_hz[0] * tau, self.epsilon);
                first_row.push((kernel_real, -kernel_imaginary_positive));
            }
            for row in 0..frequency_hz.len() {
                for column in 0..self.tau.len() {
                    let value = if row >= column {
                        first_column[row - column]
                    } else {
                        first_row[column - row]
                    };
                    real[(row, column)] = value.0;
                    imaginary[(row, column)] = value.1;
                }
            }
        } else {
            for (row, &frequency) in frequency_hz.iter().enumerate() {
                let omega = 2.0 * PI * frequency;
                for (column, &tau) in self.tau.iter().enumerate() {
                    let (kernel_real, kernel_imaginary_positive) =
                        gaussian_debye_kernel(omega * tau, self.epsilon);
                    real[(row, column)] = kernel_real;
                    imaginary[(row, column)] = -kernel_imaginary_positive;
                }
            }
        }
        Ok(BasisKernelMatrices { real, imaginary })
    }

    fn penalty(&self, order: usize, n_unpenalized: usize) -> Result<DMatrix<f64>> {
        if order != 1 {
            bail!("Gaussian basis currently supports first-order regularization only");
        }
        let n = self.tau.len();
        let mut penalty = DMatrix::zeros(n + n_unpenalized, n + n_unpenalized);
        if uses_toeplitz_tau_approximation(&self.tau) {
            let first = (0..n)
                .map(|index| {
                    gaussian_first_derivative_inner_product(
                        self.epsilon,
                        self.tau[index] / self.tau[0],
                    )
                })
                .collect::<Vec<_>>();
            for row in 0..n {
                for column in 0..n {
                    penalty[(row + n_unpenalized, column + n_unpenalized)] =
                        first[row.abs_diff(column)];
                }
            }
        } else {
            for row in 0..n {
                for column in 0..n {
                    penalty[(row + n_unpenalized, column + n_unpenalized)] =
                        gaussian_first_derivative_inner_product(
                            self.epsilon,
                            self.tau[column] / self.tau[row],
                        );
                }
            }
        }
        Ok(penalty)
    }

    fn gamma_mapping_matrix(&self) -> DMatrix<f64> {
        DMatrix::from_fn(self.tau.len(), self.tau.len(), |row, column| {
            let distance = self.epsilon * (self.tau[row] / self.tau[column]).ln();
            (-(distance * distance)).exp()
        })
    }
}

pub(super) fn evaluate_gaussian_profile(
    tau_centers: &[f64],
    coefficients: &[f64],
    epsilon: f64,
    tau_evaluation: &[f64],
) -> Result<Vec<f64>> {
    validate_tau(tau_centers)?;
    if coefficients.len() != tau_centers.len() {
        bail!(
            "Gaussian coefficient count {} does not match tau-center count {}",
            coefficients.len(),
            tau_centers.len()
        );
    }
    if !epsilon.is_finite() || epsilon <= 0.0 {
        bail!("Gaussian epsilon must be finite and strictly positive");
    }
    for (index, &tau) in tau_evaluation.iter().enumerate() {
        if !tau.is_finite() || tau <= 0.0 {
            bail!("evaluation tau[{index}] must be finite and strictly positive");
        }
    }

    Ok(tau_evaluation
        .iter()
        .map(|&tau| {
            tau_centers
                .iter()
                .zip(coefficients)
                .map(|(&center, &coefficient)| {
                    let distance = epsilon * (tau / center).ln();
                    coefficient * (-(distance * distance)).exp()
                })
                .sum()
        })
        .collect())
}

/// DRTtools' Gaussian FWHM convention:
/// `epsilon = coefficient * (2 * sqrt(ln(2))) / mean(diff(ln(tau)))`.
pub fn gaussian_epsilon_fwhm(tau: &[f64], shape_coefficient: f64) -> Result<f64> {
    validate_tau(tau)?;
    if !shape_coefficient.is_finite() || shape_coefficient <= 0.0 {
        bail!("shape coefficient must be finite and strictly positive");
    }
    let mean_log_spacing = (tau[tau.len() - 1].ln() - tau[0].ln()) / (tau.len() - 1) as f64;
    if !mean_log_spacing.is_finite() || mean_log_spacing <= 0.0 {
        bail!("cannot compute Gaussian epsilon from the tau-center spacing");
    }
    Ok(shape_coefficient * 2.0 * LN_2.sqrt() / mean_log_spacing)
}

fn validate_tau(tau: &[f64]) -> Result<()> {
    if tau.len() < 3 {
        bail!("DRT discretization requires at least three tau centers");
    }
    for (index, &value) in tau.iter().enumerate() {
        if !value.is_finite() || value <= 0.0 {
            bail!("tau[{index}] must be finite and strictly positive");
        }
        if index > 0 && value <= tau[index - 1] {
            bail!("tau grid must be strictly increasing");
        }
    }
    Ok(())
}

fn uses_drttools_toeplitz_approximation(frequency_hz: &[f64]) -> bool {
    if frequency_hz.len() < 3 {
        return false;
    }
    let spacings = frequency_hz
        .windows(2)
        .map(|pair| (1.0 / pair[1]).ln() - (1.0 / pair[0]).ln())
        .collect::<Vec<_>>();
    let mean = spacings.iter().sum::<f64>() / spacings.len() as f64;
    if mean <= 0.0 || !mean.is_finite() {
        return false;
    }
    let variance = spacings
        .iter()
        .map(|spacing| (spacing - mean).powi(2))
        .sum::<f64>()
        / (spacings.len() - 1) as f64;
    variance.sqrt() / mean < 0.01
}

fn uses_toeplitz_tau_approximation(tau: &[f64]) -> bool {
    let frequency = tau.iter().map(|tau| 1.0 / tau).collect::<Vec<_>>();
    uses_drttools_toeplitz_approximation(&frequency)
}

fn gaussian_first_derivative_inner_product(epsilon: f64, tau_ratio: f64) -> f64 {
    let a = epsilon * tau_ratio.ln();
    epsilon * (1.0 - a * a) * (-0.5 * a * a).exp() * (PI / 2.0).sqrt()
}

fn gaussian_debye_kernel(alpha: f64, epsilon: f64) -> (f64, f64) {
    let log_alpha = alpha.ln();
    let transition = -epsilon * log_alpha;
    let mut breakpoints = (-10..=10).map(f64::from).collect::<Vec<_>>();
    for point in [transition - epsilon, transition, transition + epsilon] {
        if point > -10.0 && point < 10.0 {
            breakpoints.push(point);
        }
    }
    breakpoints.sort_by(f64::total_cmp);
    breakpoints.dedup_by(|a, b| (*a - *b).abs() < 1.0e-12);

    let real = integrate_segments(&breakpoints, |t| {
        let u = log_alpha + t / epsilon;
        (-t * t).exp() * logistic_negative_twice(u) / epsilon
    });
    let imaginary = integrate_segments(&breakpoints, |t| {
        let u = (log_alpha + t / epsilon).abs();
        let exp_negative = (-u).exp();
        (-t * t).exp() * exp_negative / (1.0 + exp_negative * exp_negative) / epsilon
    });
    (real, imaginary)
}

fn logistic_negative_twice(value: f64) -> f64 {
    if value >= 0.0 {
        let exp_negative = (-2.0 * value).exp();
        exp_negative / (1.0 + exp_negative)
    } else {
        1.0 / (1.0 + (2.0 * value).exp())
    }
}

fn integrate_segments<F>(breakpoints: &[f64], function: F) -> f64
where
    F: Fn(f64) -> f64,
{
    let tolerance = 2.0e-13 / (breakpoints.len() - 1) as f64;
    breakpoints
        .windows(2)
        .map(|segment| adaptive_simpson(&function, segment[0], segment[1], tolerance, 18))
        .sum()
}

fn adaptive_simpson<F>(function: &F, left: f64, right: f64, tolerance: f64, depth: u32) -> f64
where
    F: Fn(f64) -> f64,
{
    let middle = 0.5 * (left + right);
    let f_left = function(left);
    let f_middle = function(middle);
    let f_right = function(right);
    let whole = (right - left) * (f_left + 4.0 * f_middle + f_right) / 6.0;
    adaptive_simpson_step(
        function, left, right, f_left, f_middle, f_right, whole, tolerance, depth,
    )
}

#[allow(clippy::too_many_arguments)]
fn adaptive_simpson_step<F>(
    function: &F,
    left: f64,
    right: f64,
    f_left: f64,
    f_middle: f64,
    f_right: f64,
    whole: f64,
    tolerance: f64,
    depth: u32,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let middle = 0.5 * (left + right);
    let left_middle = 0.5 * (left + middle);
    let right_middle = 0.5 * (middle + right);
    let f_left_middle = function(left_middle);
    let f_right_middle = function(right_middle);
    let left_integral = (middle - left) * (f_left + 4.0 * f_left_middle + f_middle) / 6.0;
    let right_integral = (right - middle) * (f_middle + 4.0 * f_right_middle + f_right) / 6.0;
    let refined = left_integral + right_integral;
    if depth == 0 || (refined - whole).abs() <= 15.0 * tolerance {
        return refined + (refined - whole) / 15.0;
    }
    adaptive_simpson_step(
        function,
        left,
        middle,
        f_left,
        f_left_middle,
        f_middle,
        left_integral,
        tolerance / 2.0,
        depth - 1,
    ) + adaptive_simpson_step(
        function,
        middle,
        right,
        f_middle,
        f_right_middle,
        f_right,
        right_integral,
        tolerance / 2.0,
        depth - 1,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fwhm_epsilon_matches_closed_form() {
        let tau = vec![1.0, 10.0, 100.0];
        let epsilon = gaussian_epsilon_fwhm(&tau, 0.5).unwrap();
        let expected = LN_2.sqrt() / 10.0_f64.ln();
        assert!((epsilon - expected).abs() < 1.0e-15);
    }

    #[test]
    fn gaussian_mapping_has_unit_diagonal_and_is_symmetric() {
        let basis = GaussianDiscretization::new(vec![1.0e-2, 1.0e-1, 1.0], ShapeControl::Fwhm, 0.5)
            .unwrap();
        let mapping = basis.gamma_mapping_matrix();
        for row in 0..3 {
            assert_eq!(mapping[(row, row)], 1.0);
            for column in 0..3 {
                assert!((mapping[(row, column)] - mapping[(column, row)]).abs() < 1.0e-15);
            }
        }
    }

    #[test]
    fn dense_gaussian_evaluation_matches_center_mapping() {
        let basis = GaussianDiscretization::new(vec![1.0e-2, 1.0e-1, 1.0], ShapeControl::Fwhm, 0.5)
            .unwrap();
        let coefficients = vec![0.4, 1.2, 0.7];
        let expected = basis.gamma_mapping_matrix() * DVector::from_vec(coefficients.clone());
        let actual = evaluate_gaussian_profile(
            basis.tau(),
            &coefficients,
            basis.epsilon().unwrap(),
            basis.tau(),
        )
        .unwrap();

        for (actual, expected) in actual.iter().zip(expected.iter()) {
            assert!((actual - expected).abs() < 1.0e-15);
        }
    }
}
