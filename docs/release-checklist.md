# v0.1.1 release checklist

1. Confirm a clean intended diff and the fixed MATLAB golden provenance.
2. Run, serially, formatting, Clippy, debug/release tests, docs, both release builds, `cargo package --allow-dirty`, and `git diff --check`.
3. Run `pwsh -File scripts/package_windows.ps1` and inspect the zip and `SHA256SUMS.txt`.
4. Run the packaged `eiscli.exe --version` and `--help` smoke checks.
5. Trigger `release.yml` manually and inspect its artifacts; this dry run must not create a release.
6. Commit the release state, then create and push `v0.1.1`. Only the tag event creates the GitHub Release.

Linux and macOS prebuilt archives are not part of v0.1.1.
