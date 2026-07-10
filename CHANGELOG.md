# Changelog

All notable changes to this project will be documented here. The project follows semantic versioning once a stable release process is established.

## Unreleased

### Changed

- Repaired repository ignore rules so tests, examples, and reference scripts can be committed.
- Added reproducible example-data, fixture, and golden-reference directory structure.
- Added Windows and Linux continuous integration for formatting, Clippy, tests, and release builds.
- Added baseline CLI smoke tests and documented the audited numerical limitations.
- Added a metadata-preserving EIS spectrum model with strict parsing and structured lenient-mode diagnostics.
- Added `eiscli clean`; the legacy `clean_eis` binary now delegates to the same library implementation.
- Added a deterministic fixed-size file-level batch runner with `--jobs`, failure isolation, stable output ordering, `--fail-fast`, `--resume`, and `--overwrite`.
