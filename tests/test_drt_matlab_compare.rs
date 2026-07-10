use electrochem_tools::drt::{DrtSettings, TauGridMode, solve_drt};
use electrochem_tools::drt_compare::compare_with_matlab_outputs;
use electrochem_tools::eis_io::read_eis_with_cleaning;
use std::path::Path;

#[test]
fn compares_against_exported_matlab_drttools_files() {
    let data = read_eis_with_cleaning(Path::new("examples/data/eis.z60"), true).unwrap();
    let result = solve_drt(
        &data,
        &DrtSettings {
            lambda: 1.0e-3,
            tau_min: None,
            tau_max: None,
            n_tau: 100,
            tau_grid: TauGridMode::Drttools,
            fit_inductance: false,
            regularization_order: 1,
            nonnegative: true,
            credible_intervals: false,
            solver: Default::default(),
        },
    )
    .unwrap();

    let comparison = compare_with_matlab_outputs(
        &data,
        &result,
        Some(Path::new(
            "tests/golden/drttools/eis_clean_matlab_drttools_drt_peaks.csv",
        )),
        Some(Path::new(
            "tests/golden/drttools/eis_clean_matlab_drttools_eis_regression.txt",
        )),
    )
    .unwrap();

    assert!(comparison.gamma_points_compared.unwrap() > 0);
    assert!(comparison.impedance_points_compared.unwrap() > 0);
    assert!(comparison.gamma_relative_rmse_percent.unwrap().is_finite());
    assert!(
        comparison
            .impedance_relative_rmse_percent
            .unwrap()
            .is_finite()
    );
}
