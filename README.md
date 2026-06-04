# ElectrochemTools

A collection of command-line tools for processing electrochemical workstation data, with native support for CorrTest file formats (`.cor` for CV, OCP, i-t, E-t tests and `.z60` for EIS tests).

## Tools

| Tool | Description |
|------|-------------|
| `clean_eis` | Clean EIS data by filtering out invalid data points |
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
| [`regex`](https://crates.io/crates/regex) | Pattern matching for metadata parsing |

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
