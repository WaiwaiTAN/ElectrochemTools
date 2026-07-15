use electrochem_tools::eis_io::read_eis;
use std::path::Path;

#[test]
fn reads_header_csv_and_sorts_high_to_low() {
    let data = read_eis(Path::new("examples/data/eis_simple.csv")).unwrap();
    assert_eq!(data.len(), 3);
    assert_eq!(data.frequency_hz, vec![1000.0, 100.0, 10.0]);
    assert!(data.z_imag.iter().all(|value| *value < 0.0));
}

#[test]
fn reads_clean_eis_no_header_csv() {
    let data = read_eis(Path::new("examples/data/eis_clean_eis_no_header.csv")).unwrap();
    assert_eq!(data.frequency_hz[0], 1000.0);
    assert_eq!(data.z_real[1], 2.0);
    assert_eq!(data.z_imag[1], -2.0);
}

#[test]
fn reads_raw_corrtest_z60_with_metadata() {
    let data = read_eis(Path::new("examples/data/eis.z60")).unwrap();
    assert_eq!(data.len(), 60);
    assert!(data.frequency_hz[0] > data.frequency_hz[data.len() - 1]);
    assert_eq!(data.frequency_hz[0], 100000.0);
    assert!((data.z_real[0] - 13.9207).abs() < 1.0e-4);
}

#[test]
fn reads_sanitized_bayesian_z60_fixture() {
    let data = read_eis(Path::new("tests/fixtures/bayesian_eis.z60")).unwrap();
    assert_eq!(data.len(), 70);
    assert_eq!(data.frequency_hz.first().copied(), Some(1.0e5));
    assert_eq!(data.frequency_hz.last().copied(), Some(1.0e-2));
}
