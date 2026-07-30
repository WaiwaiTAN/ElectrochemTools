# v0.1.3 release checklist

1. Confirm a clean intended diff and the fixed MATLAB golden provenance.
2. Verify `Cargo.toml`, `Cargo.lock`, `CHANGELOG.md`, and `docs/release-notes-v0.1.3.md` all identify v0.1.3.
3. Run formatting, Clippy, locked debug/release tests, `cargo package`, documentation tests, and `git diff --check`.
4. Run `pwsh -File scripts/package_windows.ps1` and inspect the standalone `.exe` files, Windows zip, and `SHA256SUMS-windows.txt`.
5. Run the packaged Windows `eiscli.exe --version`, command help, KK trimming, DRT, and ECM smoke checks.
6. Run `scripts/package_linux.sh` and `scripts/validate_linux_release.sh` in Linux/WSL; GitHub Actions repeats the native build and smoke checks.
7. Confirm the workflow validates the tag against the Cargo version and release-note path in both platform jobs.
8. Commit and push the release state to `main`, then create and push `v0.1.3`. Only the tag event creates the GitHub Release after both platform jobs succeed.

v0.1.3 publishes Windows MSVC and Linux GNU x86_64 binaries. macOS remains build-from-source only.
