use electrochem_tools::batch::{BatchOptions, BatchStatus, run_batch};
use std::fs;
use std::path::PathBuf;

fn temp_root(label: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("electrochem_tools_{label}_{}", std::process::id()));
    let _ = fs::remove_dir_all(&path);
    fs::create_dir_all(&path).unwrap();
    path
}

fn run(jobs: usize, root: PathBuf) -> Vec<String> {
    let inputs = ["a.csv", "坏 数据.csv", "c.csv"]
        .into_iter()
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    let report = run_batch(
        &inputs,
        &BatchOptions {
            jobs,
            fail_fast: false,
            overwrite: false,
            resume: false,
            out_root: root,
            output_suffix: "".to_string(),
        },
        |input, output| {
            fs::write(output.join("value.txt"), input.to_string_lossy().as_bytes())?;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(report.succeeded, 3);
    report
        .items
        .iter()
        .map(|item| fs::read_to_string(item.output_dir.join("value.txt")).unwrap())
        .collect()
}

#[test]
fn jobs_one_and_four_are_deterministic_and_support_unicode_paths() {
    assert_eq!(run(1, temp_root("serial")), run(4, temp_root("parallel")));
}

#[test]
fn a_failed_file_does_not_discard_other_results() {
    let root = temp_root("partial_failure");
    let inputs = vec![
        PathBuf::from("good-a"),
        PathBuf::from("broken"),
        PathBuf::from("good-b"),
    ];
    let report = run_batch(
        &inputs,
        &BatchOptions {
            jobs: 3,
            fail_fast: false,
            overwrite: false,
            resume: false,
            out_root: root.clone(),
            output_suffix: "".to_string(),
        },
        |input, output| {
            if input == std::path::Path::new("broken") {
                anyhow::bail!("synthetic failure");
            }
            fs::write(output.join("ok"), b"ok")?;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(report.succeeded, 2);
    assert_eq!(report.failed, 1);
    assert_eq!(report.items[1].status, BatchStatus::Failed);
    assert!(root.join("batch_summary.csv").is_file());
}
