# ElectrochemTools v0.1.3

This release adds an end-to-end Hilbert/Kramers-Kronig validation and frequency-edge trimming workflow.

`eiscli validate --out-root` can now write per-frequency directional KK residuals, a directly reusable trimmed EIS CSV, and an initial/final JSON summary. `--kk-residual-threshold` enables iterative trimming of consecutive failing points only at the high- and low-frequency edges, while `--kk-min-points` prevents over-trimming. Interior points are retained. DRT and ECM commands accept the same options and apply the retained frequency window directly before fitting.

Command help now includes copyable examples and explains the pointwise residual definition, output files, trimming behavior, and direct DRT/ECM integration.

GitHub Release assets include:

- standalone 64-bit Windows MSVC executables: `eiscli.exe`, `clean_eis.exe`, `merge_cor.exe`, and `trim_cv.exe`;
- `electrochem-tools-v0.1.3-x86_64-pc-windows-msvc.zip`, containing all Windows executables plus documentation and licenses;
- standalone Linux GNU x86_64 binaries with the `-linux-x86_64` suffix;
- `electrochem-tools-v0.1.3-x86_64-unknown-linux-gnu.tar.gz`, containing the Linux binaries plus documentation and licenses;
- `SHA256SUMS.txt`, covering every executable and archive asset.

Downloaded standalone Linux binaries may require `chmod +x <filename>` before use. The archives preserve the normal command names. macOS users should continue to build from source.

This remains a focused command-line implementation rather than a full DRTtools replacement. The KK score and trimming decisions use the documented regularized DRT cross-reconstruction proxy, not a Bayesian Hilbert Transform or Rietveld-style refinement.
