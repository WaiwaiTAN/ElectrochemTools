# DRTtools reference data

These small text exports were generated from the MATLAB DRTtools workflow and are used for diagnostic regression coverage. The local DRTtools checkout audited on 2026-07-10 was commit `034d9c4c4a4916a38a0e2f10381d931ffe1981b3`.

The current test checks readable, finite comparison metrics only. It does not compare the assembled matrices element-by-element under documented tolerances and therefore does **not** establish DRTtools parity. A later numerical-validation phase must regenerate and add `frequency.csv`, `tau.csv`, `A_re.csv`, `A_im.csv`, regularization matrices, coefficients, gamma, and reconstructed impedance before parity can be claimed.

Use `scripts/run_drttools_reference.m` with a local ignored `DRTtools/` checkout and MATLAB Optimization Toolbox to regenerate reference outputs under `target/matlab_reference/`.
