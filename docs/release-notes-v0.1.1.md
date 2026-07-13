# ElectrochemTools v0.1.1

This release expands the initial EIS toolset with DRTtools-compatible Gaussian DRT, one- and two-process RC/RQ equivalent-circuit models, and optional semi-infinite Warburg diffusion.

DRT and ECM fitting now remove positive-imaginary points by default and report the filtering count; use `--keep-positive-imag` when those points should be retained. Gaussian DRT figures evaluate the fitted RBF expansion on a dense logarithmic grid extending half a decade beyond the endpoint centers, matching the DRTtools plotting range while leaving numerical CSV outputs on the solver grid.

The archive contains 64-bit Windows MSVC builds of `eiscli`, `clean_eis`, `merge_cor`, and `trim_cv`, together with documentation and license files. Verify the archive against `SHA256SUMS.txt` before use.
