use electrochem_tools::drt::kernel::{
    assemble_combined_system, reconstruct_impedance_with_inductance,
};
use electrochem_tools::drt::regularization::piecewise_linear_penalty;
use electrochem_tools::regularization::solve_tikhonov_general_with_penalty;
use electrochem_tools::types::EisData;
use nalgebra::{DMatrix, DVector};
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    PathBuf::from("tests/golden/drttools/matrix_piecewise_linear")
}

fn matrix(path: impl AsRef<Path>) -> DMatrix<f64> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .unwrap();
    let rows = reader
        .records()
        .map(|record| {
            record
                .unwrap()
                .iter()
                .map(|value| value.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let columns = rows.first().unwrap().len();
    DMatrix::from_row_iterator(rows.len(), columns, rows.into_iter().flatten())
}

fn vector(name: &str) -> DVector<f64> {
    let values = matrix(root().join(name));
    DVector::from_iterator(values.len(), values.iter().copied())
}

fn assert_close(actual: &DMatrix<f64>, expected: &DMatrix<f64>, label: &str) {
    assert_eq!(actual.shape(), expected.shape(), "{label} shape");
    for (index, (&actual, &expected)) in actual.iter().zip(expected.iter()).enumerate() {
        let tolerance = 1.0e-10_f64.max(1.0e-8 * expected.abs());
        assert!(
            (actual - expected).abs() <= tolerance,
            "{label}[{index}]: actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
        );
    }
}

fn assert_vector_close(actual: &DVector<f64>, expected: &DVector<f64>, label: &str, relative: f64) {
    assert_eq!(actual.len(), expected.len(), "{label} length");
    for (index, (&actual, &expected)) in actual.iter().zip(expected.iter()).enumerate() {
        let tolerance = 1.0e-10_f64.max(relative * expected.abs());
        assert!(
            (actual - expected).abs() <= tolerance,
            "{label}[{index}]: actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
        );
    }
}

#[test]
fn piecewise_linear_matrices_match_matlab_drttools() {
    let frequency = vector("frequency.csv");
    let tau = vector("tau.csv");
    let observations = matrix(root().join("observations.csv"));
    let data = EisData::new(
        frequency.iter().copied().collect(),
        observations.column(0).iter().copied().collect(),
        observations.column(1).iter().copied().collect(),
    )
    .unwrap();
    let (combined, b) = assemble_combined_system(&data, tau.as_slice(), true);
    let n = data.len();
    let actual_re = combined.rows(0, n).into_owned();
    let actual_im = combined.rows(n, n).into_owned();
    assert_close(&actual_re, &matrix(root().join("A_re.csv")), "A_re");
    assert_close(&actual_im, &matrix(root().join("A_im.csv")), "A_im");

    let m1 = piecewise_linear_penalty(tau.as_slice(), 1, 2).unwrap();
    let m2 = piecewise_linear_penalty(tau.as_slice(), 2, 2).unwrap();
    assert_close(&m1, &matrix(root().join("M_1.csv")), "M_1");
    assert_close(&m2, &matrix(root().join("M_2.csv")), "M_2");

    let lambda = 1.0e-3;
    let h = (combined.transpose() * &combined + m1.scale(lambda)).scale(2.0);
    let c = (combined.transpose() * &b).scale(-2.0);
    assert_close(&h, &matrix(root().join("H.csv")), "H");
    assert_vector_close(&c, &vector("c.csv"), "c", 1.0e-8);

    let x = solve_tikhonov_general_with_penalty(&combined, &b, lambda, &m1).unwrap();
    assert_vector_close(&x, &vector("x_ridge.csv"), "x_ridge", 1.0e-7);
    assert_vector_close(
        &DVector::from_iterator(tau.len(), x.iter().skip(2).copied()),
        &vector("gamma.csv"),
        "gamma",
        1.0e-7,
    );
    let reconstructed = reconstruct_impedance_with_inductance(
        data.frequency_hz.as_slice(),
        x[1],
        x[0],
        tau.as_slice(),
        &x.iter().skip(2).copied().collect::<Vec<_>>(),
    );
    let expected_z = matrix(root().join("z_reconstructed.csv"));
    for (index, (actual, expected)) in reconstructed.iter().zip(expected_z.row_iter()).enumerate() {
        let relative_error =
            ((actual.re - expected[0]).powi(2) + (actual.im - expected[1]).powi(2)).sqrt()
                / expected[0].hypot(expected[1]).max(1.0e-300);
        assert!(
            relative_error <= 1.0e-8,
            "Z[{index}] relative error {relative_error:.3e}"
        );
    }
}
