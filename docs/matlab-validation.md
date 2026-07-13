# MATLAB numerical validation

Piecewise-linear golden data was generated on 2026-07-10 and the Gaussian Simple Run golden on 2026-07-13 with MATLAB `23.2.0.2668659 (R2023b) Update 9`, Optimization Toolbox available, and `quadprog` present. The reference checkout is Mycroft2333/DRTtools commit `034d9c4c4a4916a38a0e2f10381d931ffe1981b3`. Provenance and case configuration are embedded in every `metadata.json` and `summary.json` under `tests/golden/drttools/end_to_end`.

All cases use `examples/data/eis_cleaned.csv`, `tau = 1/frequency`, lambda `1e-3`, combined real/imaginary fitting, nonnegative gamma, and nonnegative `R_inf`. The inductance constraint is either fixed zero or nonnegative in both MATLAB and Rust.

| Case | Basis | Order | Inductance | coefficient max abs / L2 rel | gamma max abs / L2 rel | reconstructed Z max abs / L2 rel | objective rel |
| --- | --- | ---: | --- | --- | --- | --- | ---: |
| `real_order1_no_inductance` | Piecewise linear | 1 | fixed zero | `1.374e-11 / 1.253e-12` | `1.374e-11 / 1.309e-12` | `3.553e-13 / 4.139e-14` | `0` |
| `real_order1_with_inductance` | Piecewise linear | 1 | nonnegative | `1.648e-11 / 1.475e-12` | `1.648e-11 / 1.541e-12` | `9.515e-13 / 7.657e-14` | `3.409e-15` |
| `real_order2_no_inductance` | Piecewise linear | 2 | fixed zero | `6.362e-12 / 5.188e-13` | `6.362e-12 / 5.601e-13` | `1.214e-13 / 6.159e-15` | `4.267e-16` |
| `gaussian_simple_run` | Gaussian, FWHM 0.5 | 1 | fixed zero | `1.095e-12 / 1.894e-13` | `1.701e-12 / 3.598e-13` | `3.982e-13 / 2.333e-14` | `6.367e-16` |

The piecewise-linear committed tolerances are: tau `1e-12` mixed; scalar values/objective `1e-9` relative with small absolute floors; coefficients and gamma `1e-9` max absolute plus `1e-10` L2 relative; reconstructed impedance `1e-10` max absolute plus `1e-11` L2 relative. The Gaussian golden additionally compares epsilon, `A_re`, `A_im`, the analytic `M_1`, raw constrained RBF coefficients, mapped gamma, objective, and reconstructed impedance with mixed matrix tolerances and `2e-9` end-to-end bounds. Near-zero active coefficients are governed by absolute error.

Regenerate and verify locally with:

```powershell
pwsh -File scripts/regenerate_matlab_golden.ps1
```

MATLAB is not required by CI; CI reads the committed outputs. These tests validate the direct piecewise-linear formulation and the default Gaussian/combined/first-order/nonnegative/no-inductance Simple Run branch. They do not validate other RBF families, Bayesian Run/HMC, Bayesian HT, or GUI behavior.
