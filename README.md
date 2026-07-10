# ElectrochemTools

A collection of command-line tools for processing electrochemical workstation data, with native support for CorrTest file formats (`.cor` for CV, OCP, i-t, E-t tests and `.z60` for EIS tests).

## Tools

| Tool | Description |
|------|-------------|
| `clean_eis` | Clean EIS data by filtering out invalid data points |
| `eiscli` | EIS post-processing: Tikhonov DRT MVP and R(QR) equivalent-circuit fitting |
| `merge_cor` | Merge multiple `.cor` files with time-aligned timestamps |
| `trim_cv` | Trim cyclic voltammetry data to complete cycles only |

---

## Installation

### Prerequisites

- [Rust](https://www.rust-lang.org/tools/install) (toolchain supporting edition 2024)

### Build from source

```bash
git clone https://github.com/WaiwaiTAN/ElectrochemTools.git
cd ElectrochemTools
cargo build --release
```

Binaries will be placed in `target/release/`.

---

## Usage

### `clean_eis`

Filters EIS (Electrochemical Impedance Spectroscopy) data, removing rows where the imaginary impedance component (Izr) is positive — keeping only physically meaningful data points.

**Input:** Tab-separated EIS export files (9-column format, header starting with `Freq(Hz)`)

**Output:** Cleaned TSV and CSV files in the same directory as input (`{original}_cleaned.z60`, `{original}_cleaned.csv`)

```bash
# Single file
clean_eis -i eis_data.z60

# Multiple files
clean_eis -i file1.z60 file2.z60 file3.z60
```

### `eiscli`

Post-processes EIS data from cleaned CSV/TXT files or raw CorrTest `.z60` / `.txt` exports.

Subcommands:

| Command | Description |
|---------|-------------|
| `eiscli drt` | Tikhonov DRT MVP using direct Debye discretization |
| `eiscli fit-ecm` | Equivalent circuit fitting, currently `R_QR` |

Input should contain impedance columns equivalent to:

```text
frequency,Z_real,Z_imag
```

Header names such as `freq`, `frequency_hz`, `Freq(Hz)`, `Zreal`, `Z'`, `ReZ`, `Zimag`, `Z''`, `ImZ`, and CorrTest `Zr` / `Izr` are detected even when they appear after CorrTest metadata. Headerless CSV files are also supported; the first three numeric columns are interpreted as `frequency`, `Z_real`, and `Z_imag`, matching the CSV output from `clean_eis`.

Frequencies are sorted from high to low after reading. The imaginary impedance sign is preserved; use `--flip-imag` only when the input convention is known to be inverted. Use `--drop-positive-imag` to apply the same positive-imaginary filtering style as `clean_eis` while reading raw `.z60` directly.

DRT model:

```text
Z(omega) = R_inf + integral gamma(ln tau)/(1+j omega tau) dln tau
```

The default tau range is inferred as:

```text
tau_min = 1 / (2*pi*f_max) / 10
tau_max = 1 / (2*pi*f_min) * 10
```

R(QR) model:

```text
Z = Rs + 1 / (1/Rct + Q(j omega)^n)
```

Examples:

```bash
eiscli drt input.csv --lambda 1e-3 --tau-min 1e-6 --tau-max 1e3 --n-tau 100
eiscli drt eis.z60 --drop-positive-imag --auto-lambda --nonnegative --credible-intervals
eiscli drt eis.z60 --drop-positive-imag --tau-grid drttools --lambda 1e-3 --nonnegative
eiscli drt eis.z60 --drop-positive-imag --tau-grid drttools --lambda 1e-3 --nonnegative --fit-inductance
eiscli drt eis.z60 --drop-positive-imag --tau-grid drttools \
  --compare-matlab-drt eis_clean_matlab_drttools_drt_peaks.csv \
  --compare-matlab-regression eis_clean_matlab_drttools_eis_regression.txt
eiscli fit-ecm input.csv --model R_QR --out result/ --rs 0.5 --rct 20 --q 1e-3 --n 0.85
eiscli fit-ecm eis.z60 --model R_QR --auto-init --drop-positive-imag --include-correlation-matrix
```

If `--out` is omitted, output is written next to the input file using a tool-specific folder name. For example, `eis.z60` writes DRT files to `eis_drt/` and ECM fitting files to `eis_ecm/`.

Development examples:

```bash
cargo run --release --bin eiscli -- drt input.csv --lambda 1e-3 --out result/
cargo run --release --bin eiscli -- fit-ecm input.csv --model R_QR --auto-init --out result/
```

DRT outputs:

| File | Contents |
|------|----------|
| `gamma.csv` | `tau,log10_tau,gamma` |
| `drt_peaks.csv` | local DRT peak list sorted by peak height |
| `drttools_compatible_drt.csv` | MATLAB DRTtools-style `L`, `R`, then `tau,gamma(tau)` rows |
| `matlab_comparison.json` | generated when `--compare-matlab-drt` and/or `--compare-matlab-regression` is used |
| `gamma_ci.csv` | generated when `--credible-intervals` is used; linear-Gaussian 95% intervals |
| `lambda_scan.csv` | generated when `--auto-lambda` is used |
| `kk_consistency.csv` | Hilbert/Kramers-Kronig-style cross-prediction residuals |
| `kk_summary.json` | real-to-imag and imag-to-real consistency scores |
| `reconstructed_impedance.csv` | experimental impedance, fitted impedance, and residuals |
| `residual_summary.json` | lambda, tau range, nonnegative flag, RMSE, relative RMSE, and implementation note |
| `drt_gamma.svg` | DRT gamma plot |
| `nyquist_reconstruction.svg` | Nyquist plot of experimental vs reconstructed impedance |

DRT SVG plots use `log10(tau)` internally and label the x-axis at integer decades such as `10^-2`, `10^-1`, and `10^0`. Gamma plots start the y-axis at zero and use rounded tick intervals. Nyquist SVG plots use a square plotting area and equal ohm scaling on both axes, so 1 ohm horizontally has the same screen length as 1 ohm vertically.

For easier comparison with MATLAB DRTtools exports, `--tau-grid drttools` uses `tau = 1 / frequency` as the collocation grid. The default `--tau-grid logspace` keeps using a separately specified or inferred log-spaced grid.

Use `--fit-inductance` to include the DRTtools-style inductance term in the imaginary impedance model. Without this flag, inductance is fixed at zero. In `--nonnegative` mode, the bounded active-set solve uses DRTtools-compatible lower bounds for the fitted inductance, resistance, and gamma coefficients.

When MATLAB DRTtools export files are available, `--compare-matlab-drt` and `--compare-matlab-regression` produce `matlab_comparison.json` with gamma and reconstructed-impedance RMSE values. This is a diagnostic comparison, not an assertion that both tools should match exactly.

For local MATLAB reference generation without opening the DRTtools GUI, run `scripts/run_drttools_reference.m` in MATLAB batch mode from the repository root. It reads `examples/eis_cleaned.csv` by default and writes DRTtools-compatible files under `target/matlab_reference/`; the script checks for `quadprog` / Optimization Toolbox before solving.

ECM fitting outputs:

| File | Contents |
|------|----------|
| `fit_params.json` | fitted `Rs`, `Rct`, `Q`, `n`, weighted SSE, mean/reduced chi-square, RMSE, parameter standard errors, and optional correlation matrix |
| `fitted_impedance.csv` | experimental impedance, fitted impedance, and residuals |
| `nyquist_fit.svg` | Nyquist plot of experimental vs fitted impedance |

After successful `fit-ecm`, the CLI prints fit quality and parameter estimates to stdout. The parameter correlation matrix is not printed; pass `--include-correlation-matrix` to store it in `fit_params.json`.

Current limitations and TODO:

- DRT is currently a direct Debye / piecewise-linear discretization implementation, not a full reproduction of every DRTtools RBF mode.
- DRT supports a bounded active-set `--nonnegative` mode for the Tikhonov problem, but not DRTtools' MATLAB `quadprog` backend itself.
- DRT supports `--auto-lambda` scanning, local peak detection, linear-Gaussian credible intervals, and a DRT-based Hilbert/Kramers-Kronig consistency proxy.
- The credible interval output is not the original DRTtools HMC sampler; it is a deterministic Gaussian approximation around the Tikhonov solution.
- ECM fitting currently supports only `R_QR`; future models can include `Rs-(Q||R)-(Q||R)`, Warburg diffusion, and inductance.
- `modulus` and `proportional` weighting are equivalent in the current implementation because both scale complex residuals by `1 / |Z_exp|`.

Attribution: the DRT implementation in `eiscli` is derived from the algorithmic structure of the open-source MATLAB project [Mycroft2333/DRTtools](https://github.com/Mycroft2333/DRTtools), especially the real/imaginary matrix assembly, Tikhonov regularization, and EIS consistency-score workflow. This Rust CLI is an independent command-line implementation and does not include the DRTtools GUI or an exact reproduction of its RBF/QP/HMC internals.


### `merge_cor`

Merges multiple `.cor` chronoamperometry files into a single file, automatically correcting timestamps so all data aligns to a continuous time axis relative to the earliest file.

**Input:** One or more `.cor` files (tab-separated, CorrTest format with metadata header)

**Output:** A single merged `.cor` file with corrected time offsets

```bash
# Merge several sequential test files
merge_cor -i run1.cor run2.cor run3.cor -o merged.cor
```

Files are automatically sorted by test datetime before merging. The earliest file serves as the time reference; subsequent files have their time columns offset accordingly.

---

### `trim_cv`

Trims cyclic voltammetry data to retain only the last complete electrochemical cycle, discarding partial or incomplete segments.

**Input:** `.cor` CV data files containing `ExpParmas:` metadata with scan rate, frequency, and step voltages

**Output:** Trimmed file (overwrites nothing; outputs to `{original}_trimmed.{ext}`) with only complete cycle data

```bash
# Trim a single CV file
trim_cv -i cv_scan.cor

# Trim multiple files
trim_cv -i scan1.cor scan2.cor
```

**Validation:** The tool checks that the midpoint voltage of the retained cycle matches the expected step voltage (within 0.05 V tolerance), and prints warnings to stderr if a mismatch is detected.

---

## Supported File Formats

| Format | Extension | Tools |
|--------|-----------|-------|
| CorrTest CV / OCP / i-t / E-t | `.cor` | `merge_cor`, `trim_cv` |
| CorrTest EIS export | `.z60` / `.txt` | `clean_eis`, `eiscli` |


---

## Dependencies

| Crate | Purpose |
|-------|---------|
| [`clap`](https://crates.io/crates/clap) | Command-line argument parsing |
| [`anyhow`](https://crates.io/crates/anyhow) | Ergonomic error handling |
| [`calamine`](https://crates.io/crates/calamine) | Excel file reading |
| [`chrono`](https://crates.io/crates/chrono) | Date/time parsing and arithmetic |
| [`csv`](https://crates.io/crates/csv) | CSV output |
| [`nalgebra`](https://crates.io/crates/nalgebra) | Matrix solves for DRT and ECM fitting |
| [`num-complex`](https://crates.io/crates/num-complex) | Complex impedance arithmetic |
| [`regex`](https://crates.io/crates/regex) | Pattern matching for metadata parsing |
| [`serde`](https://crates.io/crates/serde) / [`serde_json`](https://crates.io/crates/serde_json) | JSON output |

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
