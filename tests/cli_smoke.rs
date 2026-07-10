use std::process::Command;

fn eiscli() -> Command {
    Command::new(env!("CARGO_BIN_EXE_eiscli"))
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
fn fit_ecm_help_succeeds() {
    assert!(
        eiscli()
            .args(["fit-ecm", "--help"])
            .status()
            .unwrap()
            .success()
    );
}
