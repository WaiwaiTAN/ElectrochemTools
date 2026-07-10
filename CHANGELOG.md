# Changelog

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
