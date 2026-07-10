$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
Push-Location $repo
try {
    matlab -batch "run('scripts/regenerate_matlab_golden.m')"
    if ($LASTEXITCODE -ne 0) { throw "MATLAB exited with code $LASTEXITCODE" }
    cargo test --locked --test test_drt_end_to_end_golden -- --nocapture
    if ($LASTEXITCODE -ne 0) { throw "Rust golden comparison failed" }
    git diff -- tests/golden/drttools/end_to_end
} finally {
    Pop-Location
}
