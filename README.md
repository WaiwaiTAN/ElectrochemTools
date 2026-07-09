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
git clone https://github.com/YOUR_USERNAME/ElectrochemTools.git
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

Post-processes EIS data from cleaned CSV/TXT files or CorrTest-style tabular exports.

Subcommands:

| Command | Description |
|---------|-------------|
| `eiscli drt` | Tikhonov DRT MVP using direct Debye discretization |
| `eiscli fit-ecm` | Equivalent circuit fitting, currently `R_QR` |

Input should contain impedance columns equivalent to:

```text
frequency,Z_real,Z_imag
```

Header names such as `freq`, `frequency_hz`, `Freq(Hz)`, `Zreal`, `Z'`, `ReZ`, `Zimag`, `Z''`, `ImZ`, and CorrTest `Zr` / `Izr` are detected. Headerless CSV files are also supported; the first three numeric columns are interpreted as `frequency`, `Z_real`, and `Z_imag`, matching the CSV output from `clean_eis`.

Frequencies are sorted from high to low after reading. The imaginary impedance sign is preserved; use `--flip-imag` only when the input convention is known to be inverted.

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
eiscli drt input.csv --lambda 1e-3 --tau-min 1e-6 --tau-max 1e3 --n-tau 100 --out result/
eiscli fit-ecm input.csv --model R_QR --out result/ --rs 0.5 --rct 20 --q 1e-3 --n 0.85
eiscli fit-ecm input.csv --model R_QR --auto-init --out result/
```

Development examples:

```bash
cargo run --release --bin eiscli -- drt input.csv --lambda 1e-3 --out result/
cargo run --release --bin eiscli -- fit-ecm input.csv --model R_QR --auto-init --out result/
```

DRT outputs:

| File | Contents |
|------|----------|
| `gamma.csv` | `tau,log10_tau,gamma` |
| `reconstructed_impedance.csv` | experimental impedance, fitted impedance, and residuals |
| `residual_summary.json` | lambda, tau range, RMSE, chi-square, and implementation note |

ECM fitting outputs:

| File | Contents |
|------|----------|
| `fit_params.json` | fitted `Rs`, `Rct`, `Q`, `n`, weighting, RMSE, chi-square, iteration count |
| `fitted_impedance.csv` | experimental impedance, fitted impedance, and residuals |

Current limitations and TODO:

- DRT is currently a direct Debye discretization MVP, not a full reproduction of the DRTtools RBF method.
- DRT gamma nonnegativity is not enforced; future work can add NNLS, bounded QP, projected gradient, or active-set solvers.
- Bayesian credible intervals, Hilbert transform / Kramers-Kronig quality scores, automatic lambda selection, and peak fitting are not implemented.
- ECM fitting currently supports only `R_QR`; future models can include `Rs-(Q||R)-(Q||R)`, Warburg diffusion, and inductance.
- `modulus` and `proportional` weighting are equivalent in the current implementation because both scale complex residuals by `1 / |Z_exp|`.
- Plot output such as `nyquist_fit.png` is a TODO.

Attribution: the DRT implementation was designed after reviewing the open-source MATLAB DRTtools algorithm structure, especially the real/imaginary matrix assembly and Tikhonov quadratic form. This Rust CLI is a smaller MVP and does not include the DRTtools GUI, RBF discretization, Bayesian analysis, or Hilbert-transform modules.


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
| CorrTest EIS export | `.z60` / `.txt` | `clean_eis` |


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
