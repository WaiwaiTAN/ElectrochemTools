# Current repository status

Audit date: 2026-07-10 (Asia/Hong_Kong)

## Repository baseline

- Branch: `main`, tracking `origin/main`.
- Starting commit: `84dcfff` (`a direct Debye / piecewise-linear discretization implementation implemented. And calculated results are compared with the matlab results.`).
- Working tree was clean before this audit.
- Toolchain used for the audit: Cargo 1.96.1 and Rust 1.96.1, `x86_64-pc-windows-msvc`.
- Local DRTtools checkout used by the reference script: commit `034d9c4c4a4916a38a0e2f10381d931ffe1981b3`. The checkout is intentionally ignored and is not distributed with this repository.

## Existing commands

| Binary | Commands or purpose |
| --- | --- |
| `eiscli` | `drt`, `fit-ecm` (`R_QR` only) |
| `clean_eis` | Legacy CorrTest EIS cleaner |
| `merge_cor` | Merge CorrTest `.cor` files |
| `trim_cv` | Retain complete CV cycles |

At the starting commit, `eiscli` did not yet have `clean` or `validate` subcommands and accepted one positional input for `drt` and `fit-ecm`.

## Baseline quality checks

The existing ignored local tests were discovered during the audit. `cargo test --all-targets` passed 7 tests, but those tests and their fixtures were excluded by `.gitignore`, so this was not reproducible from a clean clone.

| Command | Starting result |
| --- | --- |
| `cargo fmt --all -- --check` | Passed |
| `cargo clippy --all-targets --all-features -- -D warnings` | Failed with 8 diagnostics: two `too_many_arguments`, two collapsible `if` blocks, three unnecessary casts, and the non-snake-case library crate name |
| `cargo test --all-targets` | Passed locally: 7 tests; tests were ignored by Git |
| `cargo doc --no-deps` | Passed with a non-snake-case crate-name warning |

An initial attempt to run multiple Cargo quality commands concurrently caused target-directory lock contention and a transient Windows `rustc` access violation while compiling `nalgebra`. Quality commands are therefore run serially.

## Numerical status and limitations

- Deterministic direct-Debye/piecewise-linear DRT, lambda scanning, active-set non-negativity, Gaussian-linear intervals, peak detection, and a DRT cross-reconstruction score exist.
- The current MATLAB comparison test only checks that diagnostic comparison values are finite. It is **not** a matrix-level golden-reference parity test and does not establish DRTtools parity.
- The current interval implementation is a Gaussian approximation, not DRTtools HMC.
- ECM fitting supports only `R_QR`. Weighting and covariance reporting still need the statistical limitations described in the project roadmap.
- `clean_eis` and `eiscli` use different input and cleaning paths. Invalid EIS rows may currently be skipped silently in header-based input.
- There is no unified file-level batch runner, stable `run.json`, or `input_report.json` yet.

This file records the audited baseline. Later phase documents and changelog entries should distinguish implemented behavior from planned behavior.

## Repository-hygiene verification

After repairing the ignore rules and baseline diagnostics, the following commands passed serially on the audited Windows toolchain:

- `cargo fmt --all -- --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets` (10 integration tests, including 3 CLI smoke tests)
- `cargo doc --no-deps`
- `cargo build --release --all-features`
