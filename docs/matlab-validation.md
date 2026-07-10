# MATLAB numerical validation

Golden data was generated on 2026-07-10 with MATLAB `23.2.0.2668659 (R2023b) Update 9`, Optimization Toolbox available, and `quadprog` present. The reference checkout is Mycroft2333/DRTtools commit `034d9c4c4a4916a38a0e2f10381d931ffe1981b3`. Provenance and case configuration are embedded in every `metadata.json` and `summary.json` under `tests/golden/drttools/end_to_end`.

All cases use `examples/data/eis_cleaned.csv`, `tau = 1/frequency`, lambda `1e-3`, combined real/imaginary fitting, nonnegative gamma, and nonnegative `R_inf`. The inductance constraint is either fixed zero or nonnegative in both MATLAB and Rust.

| Case | Order | Inductance | coefficient max abs / L2 rel | gamma max abs / L2 rel | reconstructed Z max abs / L2 rel | objective rel |
| --- | ---: | --- | --- | --- | --- | ---: |
| `real_order1_no_inductance` | 1 | fixed zero | `1.374e-11 / 1.253e-12` | `1.374e-11 / 1.309e-12` | `3.553e-13 / 4.139e-14` | `0` |
| `real_order1_with_inductance` | 1 | nonnegative | `1.648e-11 / 1.475e-12` | `1.648e-11 / 1.541e-12` | `9.515e-13 / 7.657e-14` | `3.409e-15` |
| `real_order2_no_inductance` | 2 | fixed zero | `6.362e-12 / 5.188e-13` | `6.362e-12 / 5.601e-13` | `1.214e-13 / 6.159e-15` | `4.267e-16` |

The committed test tolerances are: tau `1e-12` mixed; scalar values/objective `1e-9` relative with small absolute floors; coefficients and gamma `1e-9` max absolute plus `1e-10` L2 relative; reconstructed impedance `1e-10` max absolute plus `1e-11` L2 relative. Near-zero active coefficients are governed by absolute error.

Regenerate and verify locally with:

```powershell
pwsh -File scripts/regenerate_matlab_golden.ps1
```

MATLAB is not required by CI; CI reads the committed outputs. These tests validate the direct piecewise-linear formulation and matched constraints, not DRTtools' RBF, HMC, Bayesian HT, or GUI behavior.
