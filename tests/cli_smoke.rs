use std::process::Command;
use std::{fs, path::PathBuf};

fn eiscli() -> Command {
    Command::new(env!("CARGO_BIN_EXE_eiscli"))
}

fn clean_eis() -> Command {
    Command::new(env!("CARGO_BIN_EXE_clean_eis"))
}

#[test]
fn top_level_help_lists_existing_commands() {
    let output = eiscli().arg("--help").output().unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("drt"));
    assert!(stdout.contains("fit-ecm"));
}

#[test]
fn drt_help_succeeds() {
    assert!(eiscli().args(["drt", "--help"]).status().unwrap().success());
}

#[test]
fn clean_help_succeeds() {
    assert!(
        eiscli()
            .args(["clean", "--help"])
            .status()
            .unwrap()
            .success()
    );
}

#[test]
fn fit_ecm_help_succeeds() {
    assert!(
        eiscli()
            .args(["fit-ecm", "--help"])
            .status()
            .unwrap()
            .success()
    );
}

#[test]
fn legacy_and_unified_clean_write_identical_data() {
    let root: PathBuf =
        std::env::temp_dir().join(format!("electrochem_tools_clean_{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let input = root.join("sample.csv");
    fs::write(
        &input,
        "frequency,z_real,z_imag\n100,1,-1\n10,2,1\n1,3,-3\n",
    )
    .unwrap();
    let unified = root.join("unified");
    let legacy = root.join("legacy");
    assert!(
        eiscli()
            .args(["clean", "-i"])
            .arg(&input)
            .arg("--out-root")
            .arg(&unified)
            .status()
            .unwrap()
            .success()
    );
    assert!(
        clean_eis()
            .arg("-i")
            .arg(&input)
            .arg("--out-root")
            .arg(&legacy)
            .status()
            .unwrap()
            .success()
    );
    assert_eq!(
        fs::read(unified.join("sample_001/cleaned.csv")).unwrap(),
        fs::read(legacy.join("sample/cleaned.csv")).unwrap()
    );
    assert!(unified.join("sample_001/input_report.json").is_file());
    assert!(legacy.join("sample/input_report.json").is_file());
}

#[test]
fn clean_batch_keeps_successful_outputs_when_one_input_fails() {
    let root =
        std::env::temp_dir().join(format!("electrochem_tools_partial_{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let good = root.join("含 空格.csv");
    let bad = root.join("bad.csv");
    fs::write(&good, "frequency,z_real,z_imag\n100,1,-1\n").unwrap();
    fs::write(&bad, "frequency,z_real,z_imag\nbroken\n").unwrap();
    let out = root.join("result");
    let status = eiscli()
        .args(["clean", "-i"])
        .arg(&good)
        .arg(&bad)
        .args(["--jobs", "2", "--out-root"])
        .arg(&out)
        .status()
        .unwrap();
    assert!(!status.success());
    assert!(out.join("sample_001/cleaned.csv").is_file());
    assert!(out.join("batch_summary.csv").is_file());
}
