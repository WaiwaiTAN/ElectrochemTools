use electrochem_tools::drt::{DrtSettings, TauGridMode, make_drttools_tau_grid, solve_drt};
use electrochem_tools::types::EisData;
use std::f64::consts::PI;

#[test]
fn drt_recovers_single_debye_peak_location() {
    let rs = 0.5;
    let r1 = 12.0;
    let tau1 = 1.0e-2;
    let freq = logspace(1.0e5, 1.0e-1, 72);
    let mut z_real = Vec::new();
    let mut z_imag = Vec::new();
    for &f in &freq {
        let wt = 2.0 * PI * f * tau1;
        z_real.push(rs + r1 / (1.0 + wt * wt));
        z_imag.push(-r1 * wt / (1.0 + wt * wt));
    }
    let data = EisData::new(freq, z_real, z_imag).unwrap();
    let result = solve_drt(
        &data,
        &DrtSettings {
            lambda: 1.0e-7,
            tau_min: Some(1.0e-5),
            tau_max: Some(1.0),
            n_tau: 120,
            tau_grid: TauGridMode::Logspace,
            basis: Default::default(),
            shape_control: Default::default(),
            shape_coefficient: 0.5,
            fit_inductance: false,
            regularization_order: 1,
            nonnegative: false,
            credible_intervals: true,
            bayesian: None,
            solver: Default::default(),
        },
    )
    .unwrap();

    assert_eq!(result.gamma.len(), 120);
    assert!(result.metrics.rmse_complex < 0.5);
    assert!(result.credible_intervals.is_some());
    assert!(result.kk.mean_score > 0.0);
    let peak_idx = result
        .gamma
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(idx, _)| idx)
        .unwrap();
    let peak_tau = result.tau[peak_idx];
    assert!((peak_tau.log10() - tau1.log10()).abs() < 0.35);
}

#[test]
fn drttools_tau_grid_uses_inverse_frequency() {
    let data = EisData::new(
        vec![100.0, 10.0, 1.0],
        vec![1.0, 2.0, 3.0],
        vec![-0.1, -0.2, -0.3],
    )
    .unwrap();
    let tau = make_drttools_tau_grid(&data).unwrap();
    assert_eq!(tau, vec![0.01, 0.1, 1.0]);
}

fn logspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    let log_start = start.log10();
    let log_end = end.log10();
    (0..n)
        .map(|idx| {
            let t = idx as f64 / (n as f64 - 1.0);
            10_f64.powf(log_start + t * (log_end - log_start))
        })
        .collect()
}
