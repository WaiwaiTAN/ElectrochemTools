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
    let output = eiscli().args(["clean", "--help"]).output().unwrap();
    assert!(output.status.success());
    assert!(
        !String::from_utf8(output.stdout)
            .unwrap()
            .contains("--resume")
    );
}

#[test]
fn fit_ecm_help_succeeds() {
    let output = eiscli().args(["fit-ecm", "--help"]).output().unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("--r2"));
    assert!(stdout.contains("--warburg-sigma"));
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
        fs::read(unified.join("sample_cleaned/cleaned.csv")).unwrap(),
        fs::read(legacy.join("sample/cleaned.csv")).unwrap()
    );
    assert!(unified.join("sample_cleaned/input_report.json").is_file());
    assert!(!unified.join("sample_cleaned/run.json").exists());
    assert!(legacy.join("sample/input_report.json").is_file());
    let report: serde_json::Value = serde_json::from_slice(
        &fs::read(unified.join("sample_cleaned/input_report.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(report["input_path"], input.display().to_string());
    assert_eq!(report["detected_input_format"], "csv");
    assert_eq!(report["mode"], "strict");
    assert_eq!(report["imag_sign"], "preserve");
    assert_eq!(report["drop_positive_imag"], true);
    assert_eq!(report["original_row_count"], 3);
    assert_eq!(report["valid_row_count"], 3);
    assert_eq!(report["removed_positive_imag_count"], 1);
    assert_eq!(report["output_point_count"], 2);
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
    assert!(out.join("含 空格_cleaned/cleaned.csv").is_file());
    assert!(out.join("batch_summary.csv").is_file());
}

#[test]
fn drt_solver_failure_returns_nonzero_without_final_result() {
    let out = std::env::temp_dir().join(format!(
        "electrochem_tools_drt_failure_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&out);
    let status = eiscli()
        .args([
            "drt",
            "-i",
            "examples/data/eis_cleaned.csv",
            "--nonnegative",
        ])
        .args([
            "--n-tau",
            "30",
            "--solver-max-iterations",
            "0",
            "--out-root",
        ])
        .arg(&out)
        .status()
        .unwrap();
    assert!(!status.success());
    assert!(!out.join("eis_drt/gamma.csv").exists());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("eis_drt/run.json")).unwrap()).unwrap();
    assert_eq!(manifest["status"], "failed");
    assert!(!out.join("eis_drt/run.json.tmp").exists());
    assert!(out.join("batch_summary.csv").is_file());
}

#[test]
fn fit_ecm_accepts_parenthesized_rc_warburg_model() {
    let out =
        std::env::temp_dir().join(format!("electrochem_tools_ecm_rc_w_{}", std::process::id()));
    let _ = fs::remove_dir_all(&out);
    assert!(
        eiscli()
            .args([
                "fit-ecm",
                "-i",
                "examples/data/eis_cleaned.csv",
                "--model",
                "R_(CR)_W",
                "--auto-init",
                "--out-root",
            ])
            .arg(&out)
            .status()
            .unwrap()
            .success()
    );
    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("eis_fit_ecm/fit_params.json")).unwrap())
            .unwrap();
    assert_eq!(summary["model"], "R_CR_W");
    assert!(summary["parameters"]["C1"].is_number());
    assert!(summary["parameters"]["sigma_w"].is_number());
}

#[test]
fn successful_drt_writes_structured_solver_report() {
    let out = std::env::temp_dir().join(format!(
        "electrochem_tools_drt_report_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&out);
    let status = eiscli()
        .args([
            "drt",
            "-i",
            "examples/data/eis_cleaned.csv",
            "--nonnegative",
        ])
        .args(["--n-tau", "30", "--out-root"])
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    let report: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("eis_drt/solver_report.json")).unwrap()).unwrap();
    assert_eq!(report["converged"], true);
    assert!(report["kkt_violation"].is_number());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("eis_drt/run.json")).unwrap()).unwrap();
    assert_eq!(manifest["status"], "success");
    assert_eq!(manifest["command"], "drt");
    assert!(manifest["input"]["sha256"].as_str().unwrap().len() == 64);
    assert!(manifest["configuration_sha256"].as_str().unwrap().len() == 64);

    let resumed = eiscli()
        .args([
            "drt",
            "-i",
            "examples/data/eis_cleaned.csv",
            "--nonnegative",
        ])
        .args(["--n-tau", "30", "--resume", "--out-root"])
        .arg(&out)
        .output()
        .unwrap();
    assert!(resumed.status.success());
    assert!(
        String::from_utf8(resumed.stdout)
            .unwrap()
            .contains("resumed")
    );

    let changed = eiscli()
        .args([
            "drt",
            "-i",
            "examples/data/eis_cleaned.csv",
            "--nonnegative",
        ])
        .args(["--n-tau", "31", "--resume", "--out-root"])
        .arg(&out)
        .output()
        .unwrap();
    assert!(!changed.status.success());
    assert!(
        String::from_utf8(changed.stderr)
            .unwrap()
            .contains("configuration differs")
    );
    assert!(
        eiscli()
            .args([
                "drt",
                "-i",
                "examples/data/eis_cleaned.csv",
                "--nonnegative",
            ])
            .args(["--n-tau", "31", "--overwrite", "--out-root"])
            .arg(&out)
            .status()
            .unwrap()
            .success()
    );
}

#[test]
fn fit_ecm_writes_and_resumes_success_manifest() {
    let out = std::env::temp_dir().join(format!(
        "electrochem_tools_ecm_manifest_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&out);
    let args = [
        "fit-ecm",
        "-i",
        "examples/data/eis_cleaned.csv",
        "--model",
        "R_QR",
        "--auto-init",
        "--out-root",
    ];
    assert!(eiscli().args(args).arg(&out).status().unwrap().success());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("eis_fit_ecm/run.json")).unwrap()).unwrap();
    assert_eq!(manifest["status"], "success");
    assert_eq!(manifest["command"], "fit-ecm");
    assert!(
        eiscli()
            .args(args)
            .arg(&out)
            .arg("--resume")
            .status()
            .unwrap()
            .success()
    );
}
