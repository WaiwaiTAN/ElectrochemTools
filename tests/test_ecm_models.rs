use electrochem_tools::ecm::{
    EcmModel, EcmModelSpec, EcmParams, EquivalentCircuitModel, ParallelBranchParams,
    ReactiveElementParams,
};
use electrochem_tools::fit::{EcmFitSettings, Weighting, fit_ecm};
use electrochem_tools::types::EisData;

#[test]
fn literature_and_compact_model_names_are_equivalent() {
    let cases = [
        ("R_(CR)", "R_CR"),
        ("R_(QR)_(CR)", "R_QR_CR"),
        ("R_(CR)_(CR)", "R_CR_CR"),
        ("R_(QR)_(QR)", "R_QR_QR"),
        ("R_(QR)_(CR)_W", "R_QR_CR_W"),
    ];
    for (input, canonical) in cases {
        let model: EcmModelSpec = input.parse().unwrap();
        assert_eq!(model.canonical_name(), canonical);
    }
}

#[test]
fn rc_and_warburg_impedance_follow_element_definitions() {
    let model: EcmModelSpec = "R_CR_W".parse().unwrap();
    let params = EcmParams {
        rs: 2.0,
        branches: vec![ParallelBranchParams {
            r: 10.0,
            reactive: ReactiveElementParams::Capacitor { c: 1.0e-3 },
        }],
        warburg_sigma: Some(3.0),
    };
    let frequency = 10.0;
    let omega = 2.0 * std::f64::consts::PI * frequency;
    let actual = EcmModel {
        spec: model,
        params,
    }
    .impedance(&[frequency])[0];
    let rc = num_complex::Complex::new(1.0 / 10.0, omega * 1.0e-3).inv();
    let warburg_scale = 3.0 / omega.sqrt();
    let expected = num_complex::Complex::new(2.0, 0.0)
        + rc
        + num_complex::Complex::new(warburg_scale, -warburg_scale);
    assert!((actual - expected).norm() < 1.0e-12);
}

#[test]
fn all_requested_topologies_fit_synthetic_impedance() {
    for name in ["R_CR", "R_QR_CR", "R_CR_CR", "R_QR_QR", "R_QR_CR_W"] {
        let spec: EcmModelSpec = name.parse().unwrap();
        let expected = parameters_for(&spec);
        let frequencies = logspace(1.0e5, 1.0e-3, 140);
        let impedance = EcmModel {
            spec: spec.clone(),
            params: expected.clone(),
        }
        .impedance(&frequencies);
        let data = EisData::new(
            frequencies,
            impedance.iter().map(|value| value.re).collect(),
            impedance.iter().map(|value| value.im).collect(),
        )
        .unwrap();
        let result = fit_ecm(
            &data,
            &EcmFitSettings {
                model: spec,
                initial: perturb(expected),
                weight: Weighting::Proportional,
                max_iter: 500,
                tol: 1.0e-12,
            },
        )
        .unwrap();
        assert!(
            result.metrics.relative_rmse_percent < 0.05,
            "{name} relative RMSE was {}%",
            result.metrics.relative_rmse_percent
        );
    }
}

fn parameters_for(spec: &EcmModelSpec) -> EcmParams {
    let branches = spec
        .branches
        .iter()
        .enumerate()
        .map(|(index, kind)| match (index, kind) {
            (0, electrochem_tools::ecm::BranchKind::Rc) => ParallelBranchParams {
                r: 8.0,
                reactive: ReactiveElementParams::Capacitor { c: 2.0e-5 },
            },
            (1, electrochem_tools::ecm::BranchKind::Rc) => ParallelBranchParams {
                r: 24.0,
                reactive: ReactiveElementParams::Capacitor { c: 8.0e-3 },
            },
            (0, electrochem_tools::ecm::BranchKind::Rq) => ParallelBranchParams {
                r: 8.0,
                reactive: ReactiveElementParams::Cpe { q: 3.0e-5, n: 0.9 },
            },
            (1, electrochem_tools::ecm::BranchKind::Rq) => ParallelBranchParams {
                r: 24.0,
                reactive: ReactiveElementParams::Cpe { q: 6.0e-3, n: 0.78 },
            },
            _ => unreachable!(),
        })
        .collect();
    EcmParams {
        rs: 0.7,
        branches,
        warburg_sigma: spec.warburg.then_some(2.5),
    }
}

fn perturb(mut params: EcmParams) -> EcmParams {
    params.rs *= 1.08;
    for branch in &mut params.branches {
        branch.r *= 0.92;
        branch.reactive = match branch.reactive {
            ReactiveElementParams::Capacitor { c } => {
                ReactiveElementParams::Capacitor { c: c * 1.12 }
            }
            ReactiveElementParams::Cpe { q, n } => ReactiveElementParams::Cpe {
                q: q * 1.12,
                n: (n - 0.03).max(0.1),
            },
        };
    }
    params.warburg_sigma = params.warburg_sigma.map(|sigma| sigma * 0.9);
    params
}

fn logspace(start: f64, end: f64, count: usize) -> Vec<f64> {
    let log_start = start.log10();
    let log_end = end.log10();
    (0..count)
        .map(|index| {
            let fraction = index as f64 / (count - 1) as f64;
            10_f64.powf(log_start + fraction * (log_end - log_start))
        })
        .collect()
}
