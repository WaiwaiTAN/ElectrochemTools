use electrochem_tools::ecm::{EquivalentCircuitModel, RqrModel, RqrParams};
use electrochem_tools::fit::{RqrFitSettings, Weighting, auto_init_rqr, fit_rqr};
use electrochem_tools::types::EisData;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};

#[test]
fn rqr_fit_recovers_synthetic_parameters() {
    let expected = RqrParams {
        rs: 0.5,
        rct: 20.0,
        q: 1.0e-3,
        n: 0.85,
    };
    let freq = logspace(1.0e5, 1.0e-1, 80);
    let model = RqrModel { params: expected };
    let z = model.impedance(&freq);
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    let noise = Normal::new(0.0, 0.002).unwrap();
    let z_real: Vec<f64> = z
        .iter()
        .map(|value| value.re * (1.0 + noise.sample(&mut rng)))
        .collect();
    let z_imag: Vec<f64> = z
        .iter()
        .map(|value| value.im * (1.0 + noise.sample(&mut rng)))
        .collect();
    let data = EisData::new(freq, z_real, z_imag).unwrap();
    let initial = auto_init_rqr(&data);
    let result = fit_rqr(
        &data,
        &RqrFitSettings {
            initial,
            weight: Weighting::Proportional,
            max_iter: 200,
            tol: 1.0e-10,
        },
    )
    .unwrap();

    assert!(relative_error(result.params.rs, expected.rs) < 0.10);
    assert!(relative_error(result.params.rct, expected.rct) < 0.15);
    assert!(relative_error(result.params.q, expected.q) < 0.30);
    assert!((result.params.n - expected.n).abs() < 0.08);
}

fn relative_error(value: f64, expected: f64) -> f64 {
    ((value - expected) / expected).abs()
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
