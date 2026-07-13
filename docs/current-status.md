# Current status

ElectrochemTools provides a strict/lenient shared EIS input layer, deterministic batch cleaning, direct-Debye/piecewise-linear Tikhonov DRT with order-1 or order-2 regularization, bounded active-set constraints, and ECM fitting for one or two RC/RQ relaxation branches with optional semi-infinite Warburg diffusion.

`eiscli clean` writes `cleaned.csv`, `cleaned.z60`, and `input_report.json`; it has no resume mode and never writes `run.json`. `eiscli drt` and `eiscli fit-ecm` write a versioned `run.json` and only resume when the command, input SHA-256, numerical configuration SHA-256, successful status, and declared outputs all match.

MATLAB R2023b/Optimization Toolbox golden validation covers the piecewise-linear matrices and three constrained end-to-end DRTtools cases at fixed upstream commit `034d9c4c4a4916a38a0e2f10381d931ffe1981b3`. See [matlab-validation.md](matlab-validation.md).

The v0.1.0 prebuilt release target is only `x86_64-pc-windows-msvc`. The core is ordinary cross-platform Rust and can be built from source elsewhere, but Linux and macOS official archives are outside this release.

This is not a full DRTtools replacement. There is no RBF DRT, HMC, Bayesian Hilbert Transform, GUI, Python binding, finite-length diffusion ECM, or inductive ECM element. Credible intervals are Gaussian linear approximations, and the KK result is a DRT cross-reconstruction proxy.
