use electrochem_tools::drt::{DrtResult, DrtSettings, TauGridMode, solve_drt};
use electrochem_tools::types::EisData;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};
use std::f64::consts::PI;

fn logspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    let (a, b) = (start.log10(), end.log10());
    (0..n)
        .map(|index| 10_f64.powf(a + index as f64 * (b - a) / (n as f64 - 1.0)))
        .collect()
}

fn spectrum(
    frequencies: Vec<f64>,
    r_inf: f64,
    inductance: f64,
    peaks: &[(f64, f64)],
    relative_noise: f64,
) -> EisData {
    let mut rng = rand::rngs::StdRng::seed_from_u64(20260710);
    let noise = Normal::new(0.0, relative_noise).unwrap();
    let mut real = Vec::with_capacity(frequencies.len());
    let mut imag = Vec::with_capacity(frequencies.len());
    for &frequency in &frequencies {
        let omega = 2.0 * PI * frequency;
        let mut re = r_inf;
        let mut im = omega * inductance;
        for &(resistance, tau) in peaks {
            let wt = omega * tau;
            re += resistance / (1.0 + wt * wt);
            im -= resistance * wt / (1.0 + wt * wt);
        }
        let scale = re.hypot(im).max(1.0);
        real.push(re + noise.sample(&mut rng) * scale);
        imag.push(im + noise.sample(&mut rng) * scale);
    }
    EisData::new(frequencies, real, imag).unwrap()
}

fn solve(data: &EisData, fit_inductance: bool, lambda: f64) -> DrtResult {
    solve_drt(
        data,
        &DrtSettings {
            lambda,
            tau_min: Some(1.0e-6),
            tau_max: Some(10.0),
            n_tau: 140,
            tau_grid: TauGridMode::Logspace,
            basis: Default::default(),
            shape_control: Default::default(),
            shape_coefficient: 0.5,
            fit_inductance,
            regularization_order: 1,
            nonnegative: true,
            credible_intervals: false,
            bayesian: None,
            solver: Default::default(),
        },
    )
    .unwrap()
}

fn has_peak_near(result: &DrtResult, expected_tau: f64, decades: f64) -> bool {
    result
        .peaks
        .iter()
        .any(|peak| (peak.tau.log10() - expected_tau.log10()).abs() <= decades)
}

#[test]
fn single_debye_recovers_resistance_peak_and_reconstruction() {
    let data = spectrum(
        logspace(1.0e5, 1.0e-1, 96),
        0.8,
        0.0,
        &[(12.0, 1.0e-2)],
        0.0,
    );
    let result = solve(&data, false, 1.0e-6);
    assert!(result.metrics.relative_rmse < 0.03);
    assert!((result.r_inf - 0.8).abs() < 0.3);
    assert!((result.polarization_resistance - 12.0).abs() / 12.0 < 0.15);
    assert!(has_peak_near(&result, 1.0e-2, 0.25));
    assert!(result.gamma.iter().all(|value| *value >= -1.0e-12));
}

#[test]
fn separated_debye_peaks_are_resolved() {
    let data = spectrum(
        logspace(1.0e6, 1.0e-2, 120),
        0.5,
        0.0,
        &[(5.0, 1.0e-4), (9.0, 1.0e-1)],
        0.0,
    );
    let result = solve(&data, false, 1.0e-6);
    assert!(result.metrics.relative_rmse < 0.04);
    assert!(has_peak_near(&result, 1.0e-4, 0.3));
    assert!(has_peak_near(&result, 1.0e-1, 0.3));
    assert!((result.polarization_resistance - 14.0).abs() / 14.0 < 0.2);
}

#[test]
fn overlapping_peaks_reconstruct_without_negative_gamma() {
    let data = spectrum(
        logspace(1.0e5, 1.0e-1, 100),
        1.0,
        0.0,
        &[(7.0, 1.0e-2), (6.0, 3.0e-2)],
        0.0,
    );
    let result = solve(&data, false, 1.0e-5);
    assert!(result.metrics.relative_rmse < 0.05);
    assert!((result.polarization_resistance - 13.0).abs() / 13.0 < 0.2);
    assert!(result.gamma.iter().all(|value| *value >= -1.0e-12));
    assert!(!result.peaks.is_empty());
}

#[test]
fn high_frequency_inductance_is_recovered_when_enabled() {
    let data = spectrum(
        logspace(1.0e5, 1.0e-1, 100),
        0.7,
        2.0e-5,
        &[(10.0, 1.0e-2)],
        0.0,
    );
    let result = solve(&data, true, 1.0e-6);
    assert!(result.metrics.relative_rmse < 0.04);
    assert!((result.inductance - 2.0e-5).abs() / 2.0e-5 < 0.2);
    assert!((result.r_inf - 0.7).abs() < 0.3);
}

#[test]
fn increasing_noise_degrades_residual_but_keeps_peak_stable() {
    let frequencies = logspace(1.0e5, 1.0e-1, 96);
    let low = solve(
        &spectrum(frequencies.clone(), 0.5, 0.0, &[(8.0, 2.0e-2)], 0.001),
        false,
        1.0e-4,
    );
    let high = solve(
        &spectrum(frequencies, 0.5, 0.0, &[(8.0, 2.0e-2)], 0.01),
        false,
        1.0e-3,
    );
    assert!(high.metrics.relative_rmse > low.metrics.relative_rmse);
    assert!(has_peak_near(&low, 2.0e-2, 0.35));
    assert!(has_peak_near(&high, 2.0e-2, 0.45));
    assert!(high.gamma.iter().all(|value| *value >= -1.0e-12));
}

#[test]
fn missing_points_and_nonuniform_grid_remain_stable() {
    let frequencies = logspace(1.0e5, 1.0e-1, 110)
        .into_iter()
        .enumerate()
        .filter_map(|(index, value)| (index % 4 != 0).then_some(value))
        .collect();
    let data = spectrum(frequencies, 0.6, 0.0, &[(11.0, 5.0e-3)], 0.002);
    let result = solve(&data, false, 1.0e-4);
    assert!(result.metrics.relative_rmse < 0.05);
    assert!(has_peak_near(&result, 5.0e-3, 0.4));
    assert!((result.r_inf - 0.6).abs() < 0.5);
}
