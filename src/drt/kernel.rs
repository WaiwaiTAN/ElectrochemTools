use crate::types::EisData;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;
use std::f64::consts::PI;

pub fn assemble_combined_system(
    data: &EisData,
    tau: &[f64],
    fit_inductance: bool,
) -> (DMatrix<f64>, DVector<f64>) {
    let weights = delta_ln_tau(tau);
    let n_points = data.len();
    let n_unpenalized = if fit_inductance { 2 } else { 1 };
    let gamma_offset = n_unpenalized;
    let r_column = if fit_inductance { 1 } else { 0 };
    let mut matrix = DMatrix::<f64>::zeros(2 * n_points, tau.len() + n_unpenalized);
    let mut observations = DVector::<f64>::zeros(2 * n_points);
    for (row, (&frequency, (&z_real, &z_imag))) in data
        .frequency_hz
        .iter()
        .zip(data.z_real.iter().zip(&data.z_imag))
        .enumerate()
    {
        let omega = 2.0 * PI * frequency;
        matrix[(row, r_column)] = 1.0;
        if fit_inductance {
            matrix[(row + n_points, 0)] = omega;
        }
        observations[row] = z_real;
        observations[row + n_points] = z_imag;
        for (column, (&tau_value, &weight)) in tau.iter().zip(&weights).enumerate() {
            let omega_tau = omega * tau_value;
            let denominator = 1.0 + omega_tau * omega_tau;
            matrix[(row, column + gamma_offset)] = weight / denominator;
            matrix[(row + n_points, column + gamma_offset)] = -weight * omega_tau / denominator;
        }
    }
    (matrix, observations)
}

pub fn assemble_real_system(data: &EisData, tau: &[f64]) -> (DMatrix<f64>, DVector<f64>) {
    let weights = delta_ln_tau(tau);
    let mut matrix = DMatrix::<f64>::zeros(data.len(), tau.len() + 1);
    let mut observations = DVector::<f64>::zeros(data.len());
    for (row, (&frequency, &z_real)) in data.frequency_hz.iter().zip(&data.z_real).enumerate() {
        let omega = 2.0 * PI * frequency;
        matrix[(row, 0)] = 1.0;
        observations[row] = z_real;
        for (column, (&tau_value, &weight)) in tau.iter().zip(&weights).enumerate() {
            let omega_tau = omega * tau_value;
            matrix[(row, column + 1)] = weight / (1.0 + omega_tau * omega_tau);
        }
    }
    (matrix, observations)
}

pub fn assemble_imag_system(data: &EisData, tau: &[f64]) -> (DMatrix<f64>, DVector<f64>) {
    let weights = delta_ln_tau(tau);
    let mut matrix = DMatrix::<f64>::zeros(data.len(), tau.len() + 1);
    let mut observations = DVector::<f64>::zeros(data.len());
    for (row, (&frequency, &z_imag)) in data.frequency_hz.iter().zip(&data.z_imag).enumerate() {
        let omega = 2.0 * PI * frequency;
        observations[row] = z_imag;
        for (column, (&tau_value, &weight)) in tau.iter().zip(&weights).enumerate() {
            let omega_tau = omega * tau_value;
            matrix[(row, column + 1)] = -weight * omega_tau / (1.0 + omega_tau * omega_tau);
        }
    }
    (matrix, observations)
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
