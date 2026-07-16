# ElectrochemTools v0.1.2

This release adds a dedicated `eiscli validate` command for lightweight Hilbert/Kramers-Kronig consistency scoring. It accepts one or more shell-expanded EIS inputs, reports directional real-to-imaginary and imaginary-to-real relative RMSE values, and presents the results in a color-aware terminal dashboard with score bars, filtering details, batch statistics, and compact paths.

DRT analysis now also includes seeded, concurrent multi-chain exact-HMC Bayesian intervals with split-R-hat and effective-sample-size diagnostics, 99% bounds, deterministic chain seeds, and resume-safe configuration. The sanitized raw `.z60` fixture and consolidated test-fixture layout support the expanded validation and Bayesian CLI coverage.

The archive contains 64-bit Windows MSVC builds of `eiscli`, `clean_eis`, `merge_cor`, and `trim_cv`, together with documentation and license files. Verify the archive against `SHA256SUMS.txt` before use.
