use electrochem_tools::eis_io::read_eis;
use electrochem_tools::types::EisData;
use std::{fs, path::Path};

fn read_temporary_csv(name: &str, contents: &str) -> EisData {
    let path = std::env::temp_dir().join(format!(
        "electrochem_tools_{name}_{}.csv",
        std::process::id()
    ));
    fs::write(&path, contents).unwrap();
    let data = read_eis(&path).unwrap();
    fs::remove_file(path).unwrap();
    data
}

#[test]
fn reads_header_csv_and_sorts_high_to_low() {
    let data = read_temporary_csv(
        "header",
        "frequency,Z_real,Z_imag\n1000,1.019999,-0.199980\n100,2,-2\n10,4.980001,-0.199980\n",
    );
    assert_eq!(data.len(), 3);
    assert_eq!(data.frequency_hz, vec![1000.0, 100.0, 10.0]);
    assert!(data.z_imag.iter().all(|value| *value < 0.0));
}

#[test]
fn reads_clean_eis_no_header_csv() {
    let data = read_temporary_csv(
        "headerless",
        "1000,1.019999,-0.199980\n100,2,-2\n10,4.980001,-0.199980\n",
    );
    assert_eq!(data.frequency_hz[0], 1000.0);
    assert_eq!(data.z_real[1], 2.0);
    assert_eq!(data.z_imag[1], -2.0);
}

#[test]
fn reads_sanitized_bayesian_z60_fixture() {
    let data = read_eis(Path::new("tests/fixtures/bayesian_eis.z60")).unwrap();
    assert_eq!(data.len(), 70);
    assert_eq!(data.frequency_hz.first().copied(), Some(1.0e5));
    assert_eq!(data.frequency_hz.last().copied(), Some(1.0e-2));
}
