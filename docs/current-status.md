# Current status

ElectrochemTools provides a strict/lenient shared EIS input layer, deterministic batch cleaning, pluggable piecewise-linear/Gaussian Tikhonov DRT discretization, bounded active-set constraints, and ECM fitting for one or two RC/RQ relaxation branches with optional semi-infinite Warburg diffusion.

`eiscli clean` and `clean_eis` both accept batches and write `<input-stem>_cleaned.csv`, `<input-stem>_cleaned.z60`, and `<input-stem>_clean_state.json` beside each input or flat below `--out-root`; cleaning has no resume mode and writes neither `batch_summary.csv` nor `run.json`. `eiscli drt` and `eiscli fit-ecm` write a versioned `run.json` and only resume when the command, input SHA-256, numerical configuration SHA-256, successful status, and declared outputs all match.

DRT and ECM fitting remove positive-imaginary points by default after optional sign flipping and report the number removed. `--keep-positive-imag` preserves them when inductive or other positive-imaginary behavior is part of the intended analysis.

MATLAB R2023b/Optimization Toolbox golden validation covers the piecewise-linear matrices, three constrained piecewise-linear end-to-end cases, and the default Gaussian Simple Run case at fixed upstream commit `034d9c4c4a4916a38a0e2f10381d931ffe1981b3`. See [matlab-validation.md](matlab-validation.md).

## DRTtools Simple Run compatibility

| Branch or option | Status |
| --- | --- |
| Combined real/imaginary fitting | Supported through the shared constrained Tikhonov backend |
| Piecewise-linear basis | Supported; regularization orders 0, 1, and 2 retain their existing behavior |
| Gaussian basis | Supported with DRTtools inverse-frequency centers and kernel integration |
| Gaussian FWHM coefficient 0.5 | Supported and MATLAB-golden validated |
| Gaussian direct shape factor | Supported, but not part of the committed default-workflow golden |
| First-order Gaussian regularization | Supported with the analytic DRTtools RBF derivative inner products |
| Nonnegative coefficients and `R_inf` | Supported by the existing bounded active-set solver |
| No-inductance Gaussian workflow | Supported and MATLAB-golden validated |
| Inductance term | Implemented in the shared backend; MATLAB goldens currently cover piecewise-linear only |
| Real-only or imaginary-only Simple Run selection | Not exposed as CLI fitting modes; used internally by the KK proxy |
| Other RBF families and Gaussian order 0/2 | Not implemented |
| Bayesian Run / HMC | Supported with concurrent deterministic multi-chain sampling, split-R-hat/ESS diagnostics, 0.5%/99.5% R-5 bounds, and nonnegative DRT coefficients; not yet MATLAB-golden validated |

For Gaussian, centers are always `tau = 1 / frequency`; piecewise-linear remains backward compatible and continues to honor `--tau-grid logspace|drttools`.

The v0.1.2 prebuilt release target is only `x86_64-pc-windows-msvc`. The core is ordinary cross-platform Rust and can be built from source elsewhere, but Linux and macOS official archives are outside this release.

This is not a full DRTtools replacement. Only Gaussian and piecewise-linear DRT bases are present. There is no Bayesian Hilbert Transform, GUI, Python binding, finite-length diffusion ECM, or inductive ECM element. Both deterministic Gaussian 95% intervals and Bayesian Run's sampled 99% bounds are available; the KK result remains a DRT cross-reconstruction proxy.
