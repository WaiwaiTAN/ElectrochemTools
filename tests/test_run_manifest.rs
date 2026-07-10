use electrochem_tools::run_manifest::{RunManifest, RunStatus, verify_resume, write_atomic};
use serde_json::json;
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);

fn fixture() -> (std::path::PathBuf, std::path::PathBuf, RunManifest) {
    let root = std::env::temp_dir().join(format!(
        "electrochem_manifest_unit_{}_{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let input = root.join("input.csv");
    fs::write(&input, "1,2,3\n").unwrap();
    let output = root.join("sample_001");
    fs::create_dir(&output).unwrap();
    fs::write(output.join("result.csv"), "ok\n").unwrap();
    let mut manifest = RunManifest::new("drt", &input, json!({"lambda": 0.001})).unwrap();
    manifest.status = RunStatus::Success;
    manifest.outputs = vec!["result.csv".into()];
    write_atomic(&output, &manifest).unwrap();
    (input, output, manifest)
}

#[test]
fn strict_resume_rejects_changed_input_missing_output_and_bad_manifest() {
    let (input, output, expected) = fixture();
    assert!(!verify_resume(&output, &expected).unwrap());

    fs::write(&input, "1,2,4\n").unwrap();
    let changed_input = RunManifest::new("drt", &input, json!({"lambda": 0.001})).unwrap();
    assert!(
        verify_resume(&output, &changed_input)
            .unwrap_err()
            .to_string()
            .contains("input SHA-256 differs")
    );
    fs::remove_file(output.join("result.csv")).unwrap();
    assert!(
        verify_resume(&output, &expected)
            .unwrap_err()
            .to_string()
            .contains("declared output is missing")
    );
    fs::write(output.join("run.json"), "not json").unwrap();
    assert!(
        verify_resume(&output, &expected)
            .unwrap_err()
            .to_string()
            .contains("invalid run.json")
    );
}

#[test]
fn strict_resume_rejects_non_success_wrong_command_schema_and_escape() {
    let (_, output, mut expected) = fixture();
    expected.status = RunStatus::Running;
    write_atomic(&output, &expected).unwrap();
    let fresh = RunManifest {
        status: RunStatus::Success,
        ..expected.clone()
    };
    assert!(
        verify_resume(&output, &fresh)
            .unwrap_err()
            .to_string()
            .contains("status is not success")
    );
    let mut existing = fresh.clone();
    existing.command = "fit-ecm".into();
    write_atomic(&output, &existing).unwrap();
    assert!(
        verify_resume(&output, &fresh)
            .unwrap_err()
            .to_string()
            .contains("command differs")
    );
    existing = fresh.clone();
    existing.schema_version = 999;
    write_atomic(&output, &existing).unwrap();
    assert!(
        verify_resume(&output, &fresh)
            .unwrap_err()
            .to_string()
            .contains("unsupported")
    );
    existing = fresh.clone();
    existing.outputs = vec!["../outside.csv".into()];
    write_atomic(&output, &existing).unwrap();
    assert!(
        verify_resume(&output, &fresh)
            .unwrap_err()
            .to_string()
            .contains("escapes")
    );
}
