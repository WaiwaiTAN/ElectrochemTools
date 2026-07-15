# Changelog

## [Unreleased]

### Added

- `eiscli drt --bayesian` now provides seeded DRTtools-style exact HMC sampling of the nonnegative DRT posterior, concurrent independent chains, split-R-hat/ESS diagnostics, 99% R-5 bounds, Bayesian CSV/JSON/SVG outputs, and resume-safe sampler configuration.
- Added a metadata-sanitized real `.z60` fixture for Bayesian CLI and input-reader tests.

### Changed

- Consolidated sample inputs under `tests/fixtures` and removed the redundant generated `examples` tree.

## [0.1.1] - 2026-07-13

### Added

- ECM models `R_CR`, `R_QR_CR`, `R_CR_CR`, and `R_QR_QR`, with an optional series semi-infinite Warburg element selected by the `_W` suffix.
- Numbered branch initialization options and dynamic parameter/uncertainty labels in `fit_params.json`; the original `R_QR` model and its CLI option aliases remain supported.

### Changed

- DRT and ECM fitting now drop positive-imaginary points by default, report the filtering count to the console, and provide `--keep-positive-imag` to retain them.
- Gaussian DRT SVG plots now evaluate the fitted RBF expansion on a 10x dense logarithmic grid extending half a decade beyond both endpoint centers, while preserving solver-grid CSV outputs.
- SVG plots now store LaTeX-ready axis labels with `siunitx` units and semantic text IDs directly in the single viewable SVG output.

### Fixed

- `eiscli clean` and `clean_eis` now share batch cleaning and write `<input>_cleaned.csv`, `<input>_cleaned.z60`, and `<input>_clean_state.json` beside each source (or flat below `--out-root`) instead of creating per-input result directories.

## [0.1.0] - 2026-07-10

### Added

- Unified strict/lenient EIS input, cleaning diagnostics, stable batch output, direct-Debye/piecewise-linear DRT, and `R_QR` fitting.
- Order-1/order-2 regularization, bounded active-set constraints, Gaussian interval estimates, and KK cross-reconstruction diagnostics.
- Matrix-level and three constrained end-to-end MATLAB R2023b/DRTtools golden cases with numerical error thresholds and provenance.
- Strict SHA-256-based `run.json` resume for DRT and ECM calculations.
- Reproducible Windows MSVC packaging, checksums, and tag/manual GitHub Actions workflow.

### Changed

- `eiscli clean` no longer accepts `--resume`; it writes `input_report.json` but no computation manifest.
- Batch summaries distinguish successful, failed, resumed, and not-processed inputs.

### Fixed

- Resume no longer treats an existing output directory as proof of a completed matching calculation.
- Clean reports now expose source format, cleaning policy, row counts, skip reasons, filtering count, and output files.

### Known limitations

- This is not a full DRTtools replacement: no RBF, HMC, Bayesian HT, or MATLAB GUI.
- ECM fitting supports only `R_QR`; intervals are Gaussian approximations and KK is a cross-reconstruction proxy.
- Official v0.1.0 binaries are Windows x86-64 MSVC only.
