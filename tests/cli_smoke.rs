use std::process::Command;
use std::{fs, path::PathBuf};

const BAYESIAN_FIXTURE: &str = "tests/fixtures/bayesian_eis.z60";

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
    let output = eiscli().args(["drt", "--help"]).output().unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("--keep-positive-imag"));
    assert!(stdout.contains("--basis"));
    assert!(stdout.contains("--shape-control"));
    assert!(stdout.contains("--shape-coefficient"));
    assert!(stdout.contains("--bayesian"));
    assert!(stdout.contains("--bayesian-run"));
    assert!(stdout.contains("--bayesian-samples"));
    assert!(stdout.contains("--bayesian-burn-in"));
    assert!(stdout.contains("--bayesian-seed"));
    assert!(stdout.contains("--bayesian-chains"));
    assert!(!stdout.contains("--drop-positive-imag"));
}

#[test]
fn bayesian_drt_writes_seeded_hmc_outputs_and_manifest_configuration() {
    let out = std::env::temp_dir().join(format!(
        "electrochem_tools_bayesian_drt_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&out);
    let output = eiscli()
        .args([
            "drt",
            "-i",
            BAYESIAN_FIXTURE,
            "--n-tau",
            "3",
            "--bayesian",
            "--bayesian-samples",
            "1000",
            "--bayesian-burn-in",
            "100",
            "--bayesian-seed",
            "42",
            "--bayesian-chains",
            "4",
            "--out-root",
        ])
        .arg(&out)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let result = out.join("bayesian_eis_drt");
    let csv = fs::read_to_string(result.join("gamma_bayesian.csv")).unwrap();
    let mut rows = csv.lines();
    assert_eq!(
        rows.next().unwrap(),
        "tau,log10_tau,gamma_map,gamma_mean,gamma_lower_99,gamma_upper_99"
    );
    assert_eq!(rows.count(), 3);
    assert!(result.join("bayesian_summary.json").is_file());
    assert!(result.join("bayesian_diagnostics.csv").is_file());
    assert!(result.join("drt_gamma_bayesian.svg").is_file());

    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(result.join("bayesian_summary.json")).unwrap()).unwrap();
    assert_eq!(summary["method"], "exact_constrained_gaussian_hmc");
    assert_eq!(summary["chains"], 4);
    assert_eq!(summary["samples_per_chain"], 1000);
    assert_eq!(summary["total_draws"], 4000);
    assert_eq!(summary["burn_in_per_chain"], 100);
    assert_eq!(summary["retained_draws_per_chain"], 900);
    assert_eq!(summary["retained_draws"], 3600);
    assert_eq!(summary["seed"], 42);
    assert_eq!(summary["chain_seeds"].as_array().unwrap().len(), 4);
    assert_eq!(summary["credible_level"], 0.99);
    assert!(summary["max_split_r_hat"].is_number());
    assert!(summary["min_effective_sample_size"].is_number());

    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(result.join("run.json")).unwrap()).unwrap();
    assert_eq!(manifest["configuration"]["bayesian"]["enabled"], true);
    assert_eq!(manifest["configuration"]["bayesian"]["chains"], 4);
    assert_eq!(manifest["configuration"]["nonnegative"], true);
    assert!(
        manifest["outputs"]
            .as_array()
            .unwrap()
            .iter()
            .any(|path| path == "gamma_bayesian.csv")
    );
}

#[test]
fn bayesian_drt_rejects_invalid_or_misplaced_sampler_options() {
    let missing_mode = eiscli()
        .args(["drt", "-i", BAYESIAN_FIXTURE, "--bayesian-seed", "42"])
        .output()
        .unwrap();
    assert!(!missing_mode.status.success());

    let too_few = eiscli()
        .args([
            "drt",
            "-i",
            BAYESIAN_FIXTURE,
            "--bayesian",
            "--bayesian-samples",
            "999",
        ])
        .output()
        .unwrap();
    assert!(!too_few.status.success());
    assert!(
        String::from_utf8_lossy(&too_few.stderr)
            .contains("requires at least 1000 samples per chain")
    );

    let conflicting_intervals = eiscli()
        .args([
            "drt",
            "-i",
            BAYESIAN_FIXTURE,
            "--bayesian",
            "--credible-intervals",
        ])
        .output()
        .unwrap();
    assert!(!conflicting_intervals.status.success());
}

#[test]
fn gaussian_drt_cli_uses_drttools_centers_and_records_shape_settings() {
    let out = std::env::temp_dir().join(format!(
        "electrochem_tools_gaussian_drt_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&out);
    let status = eiscli()
        .args([
            "drt",
            "-i",
            BAYESIAN_FIXTURE,
            "--basis",
            "gaussian",
            "--shape-control",
            "fwhm",
            "--shape-coefficient",
            "0.5",
            "--lambda",
            "1e-3",
            "--regularization-order",
            "1",
            "--nonnegative",
            "--out-root",
        ])
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    let summary: serde_json::Value = serde_json::from_slice(
        &fs::read(out.join("bayesian_eis_drt/residual_summary.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(summary["basis"], "gaussian");
    assert_eq!(summary["tau_grid"], "drttools");
    assert_eq!(summary["shape_control"], "fwhm");
    assert_eq!(summary["shape_coefficient"], 0.5);
    assert!(summary["epsilon"].as_f64().unwrap() > 0.0);
    assert_eq!(summary["n_tau"], summary["n_points"]);

    let gamma_rows = fs::read_to_string(out.join("bayesian_eis_drt/gamma.csv"))
        .unwrap()
        .lines()
        .count()
        - 1;
    let svg = fs::read_to_string(out.join("bayesian_eis_drt/drt_gamma.svg")).unwrap();
    assert!(svg.contains(r#"id="title""#));
    assert!(svg.contains(r#"id="xlabel""#));
    assert!(svg.contains(r"$\tau$ / \unit{\second}"));
    assert!(svg.contains(r#"id="ylabel""#));
    assert!(svg.contains(r"$\gamma(\ln\tau)$ / \unit{\ohm}"));
    assert!(svg.contains(r#"id="legend-gamma""#));
    assert!(svg.contains("$10^{-1}$"));
    assert!(!svg.contains(">10^-1<"));
    assert!(!out.join("bayesian_eis_drt/drt_gamma_latex.svg").exists());
    assert!(!out.join("bayesian_eis_drt/drt_gamma_plain.svg").exists());
    let plotted_points = svg
        .split_once("<polyline")
        .unwrap()
        .1
        .split_once("points=\"")
        .unwrap()
        .1
        .split_once('"')
        .unwrap()
        .0
        .split_whitespace()
        .count();
    assert_eq!(plotted_points, 10 * gamma_rows);
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
    assert!(stdout.contains("--keep-positive-imag"));
    assert!(!stdout.contains("--drop-positive-imag"));
}

#[test]
fn analysis_filters_positive_imag_by_default_and_reports_count() {
    let root = std::env::temp_dir().join(format!(
        "electrochem_tools_positive_imag_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let input = root.join("spectrum.csv");
    fs::write(
        &input,
        concat!(
            "frequency,z_real,z_imag\n",
            "10000,1.0,-0.01\n",
            "1000,1.1,-0.1\n",
            "100,2.0,-1.0\n",
            "10,4.0,-2.0\n",
            "1,5.0,-0.5\n",
            "0.1,5.5,0.2\n",
        ),
    )
    .unwrap();

    let fit_root = root.join("fit");
    let fit = eiscli()
        .args(["fit-ecm", "-i"])
        .arg(&input)
        .args(["--model", "R_CR", "--auto-init", "--out-root"])
        .arg(&fit_root)
        .output()
        .unwrap();
    assert!(
        fit.status.success(),
        "{}",
        String::from_utf8_lossy(&fit.stderr)
    );
    assert!(
        String::from_utf8(fit.stdout)
            .unwrap()
            .contains("positive-imaginary filter (ECM fitting): removed 1 of 6 point(s)")
    );
    let fit_summary: serde_json::Value = serde_json::from_slice(
        &fs::read(fit_root.join("spectrum_fit_ecm/fit_params.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(fit_summary["n_points"], 5);

    let drt_root = root.join("drt");
    let drt = eiscli()
        .args(["drt", "-i"])
        .arg(&input)
        .args(["--n-tau", "8", "--out-root"])
        .arg(&drt_root)
        .output()
        .unwrap();
    assert!(
        drt.status.success(),
        "{}",
        String::from_utf8_lossy(&drt.stderr)
    );
    assert!(
        String::from_utf8(drt.stdout)
            .unwrap()
            .contains("positive-imaginary filter (DRT): removed 1 of 6 point(s)")
    );

    let keep_root = root.join("keep");
    let keep = eiscli()
        .args(["fit-ecm", "-i"])
        .arg(&input)
        .args([
            "--model",
            "R_CR",
            "--auto-init",
            "--keep-positive-imag",
            "--out-root",
        ])
        .arg(&keep_root)
        .output()
        .unwrap();
    assert!(
        keep.status.success(),
        "{}",
        String::from_utf8_lossy(&keep.stderr)
    );
    assert!(
        !String::from_utf8(keep.stdout)
            .unwrap()
            .contains("positive-imaginary filter")
    );
    let keep_summary: serde_json::Value = serde_json::from_slice(
        &fs::read(keep_root.join("spectrum_fit_ecm/fit_params.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(keep_summary["n_points"], 6);
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
        fs::read(unified.join("sample_cleaned.csv")).unwrap(),
        fs::read(legacy.join("sample_cleaned.csv")).unwrap()
    );
    assert!(unified.join("sample_cleaned.z60").is_file());
    assert!(unified.join("sample_clean_state.json").is_file());
    assert!(legacy.join("sample_cleaned.z60").is_file());
    assert!(legacy.join("sample_clean_state.json").is_file());
    assert!(!unified.join("sample_cleaned").exists());
    assert!(!unified.join("batch_summary.csv").exists());
    assert!(!legacy.join("batch_summary.csv").exists());
    let report: serde_json::Value =
        serde_json::from_slice(&fs::read(unified.join("sample_clean_state.json")).unwrap())
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
fn eiscli_clean_defaults_to_three_files_beside_the_input() {
    let root = std::env::temp_dir().join(format!(
        "electrochem_tools_sibling_clean_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let input = root.join("spectrum.csv");
    fs::write(&input, "frequency,z_real,z_imag\n100,1,-1\n").unwrap();

    assert!(
        eiscli()
            .args(["clean", "-i"])
            .arg(&input)
            .status()
            .unwrap()
            .success()
    );
    assert!(root.join("spectrum_cleaned.csv").is_file());
    assert!(root.join("spectrum_cleaned.z60").is_file());
    assert!(root.join("spectrum_clean_state.json").is_file());
    assert!(!root.join("spectrum_cleaned").exists());
    assert!(!root.join("batch_summary.csv").exists());
}

#[test]
fn clean_eis_cleans_multiple_inputs_in_one_batch() {
    let root = std::env::temp_dir().join(format!(
        "electrochem_tools_legacy_batch_clean_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let first = root.join("first.csv");
    let second = root.join("second.csv");
    fs::write(&first, "frequency,z_real,z_imag\n100,1,-1\n").unwrap();
    fs::write(&second, "frequency,z_real,z_imag\n10,2,-2\n").unwrap();

    assert!(
        clean_eis()
            .arg("-i")
            .arg(&first)
            .arg(&second)
            .args(["--jobs", "2"])
            .status()
            .unwrap()
            .success()
    );
    for stem in ["first", "second"] {
        assert!(root.join(format!("{stem}_cleaned.csv")).is_file());
        assert!(root.join(format!("{stem}_cleaned.z60")).is_file());
        assert!(root.join(format!("{stem}_clean_state.json")).is_file());
    }
    assert!(!root.join("batch_summary.csv").exists());
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
    assert!(out.join("含 空格_cleaned.csv").is_file());
    assert!(out.join("含 空格_cleaned.z60").is_file());
    assert!(out.join("含 空格_clean_state.json").is_file());
    assert!(!out.join("batch_summary.csv").exists());
}

#[test]
fn drt_solver_failure_returns_nonzero_without_final_result() {
    let out = std::env::temp_dir().join(format!(
        "electrochem_tools_drt_failure_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&out);
    let status = eiscli()
        .args(["drt", "-i", BAYESIAN_FIXTURE, "--nonnegative"])
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
    assert!(!out.join("bayesian_eis_drt/gamma.csv").exists());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("bayesian_eis_drt/run.json")).unwrap()).unwrap();
    assert_eq!(manifest["status"], "failed");
    assert!(!out.join("bayesian_eis_drt/run.json.tmp").exists());
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
                BAYESIAN_FIXTURE,
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
    let summary: serde_json::Value = serde_json::from_slice(
        &fs::read(out.join("bayesian_eis_fit_ecm/fit_params.json")).unwrap(),
    )
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
        .args(["drt", "-i", BAYESIAN_FIXTURE, "--nonnegative"])
        .args(["--n-tau", "30", "--out-root"])
        .arg(&out)
        .status()
        .unwrap();
    assert!(status.success());
    let report: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("bayesian_eis_drt/solver_report.json")).unwrap())
            .unwrap();
    assert_eq!(report["converged"], true);
    assert!(report["kkt_violation"].is_number());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("bayesian_eis_drt/run.json")).unwrap()).unwrap();
    assert_eq!(manifest["status"], "success");
    assert_eq!(manifest["command"], "drt");
    assert!(manifest["input"]["sha256"].as_str().unwrap().len() == 64);
    assert!(manifest["configuration_sha256"].as_str().unwrap().len() == 64);

    let resumed = eiscli()
        .args(["drt", "-i", BAYESIAN_FIXTURE, "--nonnegative"])
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
        .args(["drt", "-i", BAYESIAN_FIXTURE, "--nonnegative"])
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
            .args(["drt", "-i", BAYESIAN_FIXTURE, "--nonnegative",])
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
        BAYESIAN_FIXTURE,
        "--model",
        "R_QR",
        "--auto-init",
        "--out-root",
    ];
    assert!(eiscli().args(args).arg(&out).status().unwrap().success());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("bayesian_eis_fit_ecm/run.json")).unwrap())
            .unwrap();
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
