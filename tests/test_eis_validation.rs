use electrochem_tools::eis::{
    ImagSignPolicy, ImaginaryConvention, ReadOptions, apply_imag_sign, read_spectrum,
};
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_ID: AtomicU64 = AtomicU64::new(0);

fn fixture(contents: &str) -> PathBuf {
    let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("electrochem_tools_{}_{}", std::process::id(), id));
    fs::create_dir_all(&dir).unwrap();
    let path = dir.join("input.csv");
    fs::write(&path, contents).unwrap();
    path
}

#[test]
fn strict_reader_rejects_non_finite_and_zero_frequency() {
    for row in ["NaN,1,-1", "100,Inf,-1", "0,1,-1"] {
        let path = fixture(&format!("frequency,z_real,z_imag\n{row}\n"));
        assert!(read_spectrum(&path, &ReadOptions::default()).is_err());
    }
}

#[test]
fn strict_reader_rejects_bad_rows_and_duplicate_frequency() {
    let bad = fixture("frequency,z_real,z_imag\n100,1,-1\nbroken,row\n");
    assert!(read_spectrum(&bad, &ReadOptions::default()).is_err());

    let duplicate = fixture("frequency,z_real,z_imag\n100,1,-1\n100,2,-2\n");
    assert!(read_spectrum(&duplicate, &ReadOptions::default()).is_err());
}

#[test]
fn lenient_reader_reports_skipped_rows_and_preserves_sign() {
    let path = fixture("frequency,z_real,z_imag\n100,1,1\nbroken,row\n10,2,2\n0,3,3\n");
    let parsed = read_spectrum(&path, &ReadOptions { lenient: true }).unwrap();
    assert_eq!(parsed.spectrum.points.len(), 2);
    assert_eq!(parsed.report.rows_read, 2);
    assert_eq!(parsed.report.rows_skipped, 2);
    assert_eq!(
        parsed.spectrum.metadata.imaginary_convention,
        ImaginaryConvention::PositiveCapacitive
    );
    assert!(
        parsed
            .spectrum
            .points
            .iter()
            .all(|point| point.z_imag_ohm > 0.0)
    );
}

#[test]
fn strict_reader_rejects_ambiguous_headers() {
    let path = fixture("frequency,freq,z_real,z_imag\n100,100,1,-1\n");
    assert!(read_spectrum(&path, &ReadOptions::default()).is_err());
}

#[test]
fn imaginary_sign_changes_only_for_explicit_policy() {
    let path = fixture("frequency,z_real,z_imag\n100,1,1\n10,2,2\n");
    let mut spectrum = read_spectrum(&path, &ReadOptions::default())
        .unwrap()
        .spectrum;
    apply_imag_sign(&mut spectrum, ImagSignPolicy::Preserve);
    assert!(spectrum.points.iter().all(|point| point.z_imag_ohm > 0.0));
    apply_imag_sign(&mut spectrum, ImagSignPolicy::NegativeCapacitive);
    assert!(spectrum.points.iter().all(|point| point.z_imag_ohm < 0.0));
}
