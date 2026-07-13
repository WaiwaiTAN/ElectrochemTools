use crate::types::EisData;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct BasisKernelMatrices {
    pub real: DMatrix<f64>,
    pub imaginary: DMatrix<f64>,
}

pub fn piecewise_linear_kernel_matrices(frequency_hz: &[f64], tau: &[f64]) -> BasisKernelMatrices {
    let weights = delta_ln_tau(tau);
    let mut real = DMatrix::<f64>::zeros(frequency_hz.len(), tau.len());
    let mut imaginary = DMatrix::<f64>::zeros(frequency_hz.len(), tau.len());
    for (row, &frequency) in frequency_hz.iter().enumerate() {
        let omega = 2.0 * PI * frequency;
        for (column, (&tau_value, &weight)) in tau.iter().zip(&weights).enumerate() {
            let omega_tau = omega * tau_value;
            let denominator = 1.0 + omega_tau * omega_tau;
            real[(row, column)] = weight / denominator;
            imaginary[(row, column)] = -weight * omega_tau / denominator;
        }
    }
    BasisKernelMatrices { real, imaginary }
}

pub fn assemble_combined_from_kernels(
    data: &EisData,
    kernels: &BasisKernelMatrices,
    fit_inductance: bool,
) -> (DMatrix<f64>, DVector<f64>) {
    let n_points = data.len();
    let n_unpenalized = if fit_inductance { 2 } else { 1 };
    let gamma_offset = n_unpenalized;
    let r_column = if fit_inductance { 1 } else { 0 };
    let mut matrix = DMatrix::<f64>::zeros(2 * n_points, kernels.real.ncols() + n_unpenalized);
    let mut observations = DVector::<f64>::zeros(2 * n_points);
    for (row, (&frequency, (&z_real, &z_imag))) in data
        .frequency_hz
        .iter()
        .zip(data.z_real.iter().zip(&data.z_imag))
        .enumerate()
    {
        matrix[(row, r_column)] = 1.0;
        if fit_inductance {
            matrix[(row + n_points, 0)] = 2.0 * PI * frequency;
        }
        observations[row] = z_real;
        observations[row + n_points] = z_imag;
        for column in 0..kernels.real.ncols() {
            matrix[(row, column + gamma_offset)] = kernels.real[(row, column)];
            matrix[(row + n_points, column + gamma_offset)] = kernels.imaginary[(row, column)];
        }
    }
    (matrix, observations)
}

pub fn assemble_real_from_kernels(
    data: &EisData,
    kernels: &BasisKernelMatrices,
) -> (DMatrix<f64>, DVector<f64>) {
    let mut matrix = DMatrix::<f64>::zeros(data.len(), kernels.real.ncols() + 1);
    let observations = DVector::from_iterator(data.len(), data.z_real.iter().copied());
    for row in 0..data.len() {
        matrix[(row, 0)] = 1.0;
        for column in 0..kernels.real.ncols() {
            matrix[(row, column + 1)] = kernels.real[(row, column)];
        }
    }
    (matrix, observations)
}

pub fn assemble_imag_from_kernels(
    data: &EisData,
    kernels: &BasisKernelMatrices,
) -> (DMatrix<f64>, DVector<f64>) {
    let mut matrix = DMatrix::<f64>::zeros(data.len(), kernels.imaginary.ncols() + 1);
    let observations = DVector::from_iterator(data.len(), data.z_imag.iter().copied());
    for row in 0..data.len() {
        for column in 0..kernels.imaginary.ncols() {
            matrix[(row, column + 1)] = kernels.imaginary[(row, column)];
        }
    }
    (matrix, observations)
}

pub fn reconstruct_from_kernels(
    frequency_hz: &[f64],
    r_inf: f64,
    inductance: f64,
    kernels: &BasisKernelMatrices,
    coefficients: &[f64],
) -> Vec<Complex<f64>> {
    assert_eq!(frequency_hz.len(), kernels.real.nrows());
    assert_eq!(kernels.real.shape(), kernels.imaginary.shape());
    assert_eq!(coefficients.len(), kernels.real.ncols());
    (0..frequency_hz.len())
        .map(|row| {
            let mut value = Complex::new(r_inf, 2.0 * PI * frequency_hz[row] * inductance);
            for (column, &coefficient) in coefficients.iter().enumerate() {
                value.re += kernels.real[(row, column)] * coefficient;
                value.im += kernels.imaginary[(row, column)] * coefficient;
            }
            value
        })
        .collect()
}

pub fn assemble_combined_system(
    data: &EisData,
    tau: &[f64],
    fit_inductance: bool,
) -> (DMatrix<f64>, DVector<f64>) {
    let kernels = piecewise_linear_kernel_matrices(&data.frequency_hz, tau);
    assemble_combined_from_kernels(data, &kernels, fit_inductance)
}

pub fn assemble_real_system(data: &EisData, tau: &[f64]) -> (DMatrix<f64>, DVector<f64>) {
    let kernels = piecewise_linear_kernel_matrices(&data.frequency_hz, tau);
    assemble_real_from_kernels(data, &kernels)
}

pub fn assemble_imag_system(data: &EisData, tau: &[f64]) -> (DMatrix<f64>, DVector<f64>) {
    let kernels = piecewise_linear_kernel_matrices(&data.frequency_hz, tau);
    assemble_imag_from_kernels(data, &kernels)
}

pub fn reconstruct_impedance(
    frequency_hz: &[f64],
    r_inf: f64,
    tau: &[f64],
    gamma: &[f64],
) -> Vec<Complex<f64>> {
    reconstruct_impedance_with_inductance(frequency_hz, r_inf, 0.0, tau, gamma)
}

pub fn reconstruct_impedance_with_inductance(
    frequency_hz: &[f64],
    r_inf: f64,
    inductance: f64,
    tau: &[f64],
    gamma: &[f64],
) -> Vec<Complex<f64>> {
    let weights = delta_ln_tau(tau);
    frequency_hz
        .iter()
        .map(|&frequency| {
            let omega = 2.0 * PI * frequency;
            tau.iter().zip(gamma).zip(&weights).fold(
                Complex::new(r_inf, omega * inductance),
                |impedance, ((&tau_value, &gamma_value), &weight)| {
                    let omega_tau = omega * tau_value;
                    let denominator = 1.0 + omega_tau * omega_tau;
                    impedance
                        + Complex::new(
                            gamma_value * weight / denominator,
                            -gamma_value * weight * omega_tau / denominator,
                        )
                },
            )
        })
        .collect()
}

pub fn delta_ln_tau(tau: &[f64]) -> Vec<f64> {
    if tau.len() < 2 {
        return vec![1.0; tau.len()];
    }
    let logarithms: Vec<f64> = tau.iter().map(|value| value.ln()).collect();
    (0..tau.len())
        .map(|index| {
            if index == 0 {
                0.5 * (logarithms[1] - logarithms[0])
            } else if index + 1 == tau.len() {
                0.5 * (logarithms[index] - logarithms[index - 1])
            } else {
                0.5 * (logarithms[index + 1] - logarithms[index - 1])
            }
        })
        .collect()
}
