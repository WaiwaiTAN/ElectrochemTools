param([string]$Version)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$manifest = Join-Path $repo 'Cargo.toml'
if (-not $Version) {
    $match = Select-String -Path $manifest -Pattern '^version\s*=\s*"([^"]+)"' | Select-Object -First 1
    if (-not $match) { throw 'Could not read package version from Cargo.toml' }
    $Version = $match.Matches[0].Groups[1].Value
}
$target = 'x86_64-pc-windows-msvc'
$name = "electrochem-tools-v$Version-$target"
$dist = Join-Path $repo 'target\dist'
$stage = Join-Path $repo "target\package-stage\$name"
$zip = Join-Path $dist "$name.zip"
$sums = Join-Path $dist 'SHA256SUMS.txt'

Push-Location $repo
try {
    $built = $false
    for ($attempt = 1; $attempt -le 3; $attempt++) {
        cargo build --locked --release --all-features --target $target --jobs 1
        if ($LASTEXITCODE -eq 0) { $built = $true; break }
        Write-Warning "release build attempt $attempt failed; retrying transient compiler failure"
    }
    if (-not $built) { throw 'cargo build failed after 3 attempts' }
    $resolvedRepo = (Resolve-Path $repo).Path
    foreach ($path in @($stage, $zip, $sums)) {
        $full = [System.IO.Path]::GetFullPath($path)
        if (-not $full.StartsWith($resolvedRepo, [System.StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to remove path outside repository: $full"
        }
        if (Test-Path -LiteralPath $full) { Remove-Item -LiteralPath $full -Recurse -Force }
    }
    New-Item -ItemType Directory -Force -Path $stage, $dist | Out-Null
    foreach ($binary in @('eiscli.exe', 'clean_eis.exe', 'merge_cor.exe', 'trim_cv.exe')) {
        $source = Join-Path $repo "target\$target\release\$binary"
        if (-not (Test-Path -LiteralPath $source)) { throw "Missing binary: $source" }
        Copy-Item -LiteralPath $source -Destination $stage
    }
    foreach ($document in @('README.md', 'LICENSE', 'THIRD_PARTY_NOTICES.md', 'CHANGELOG.md')) {
        Copy-Item -LiteralPath (Join-Path $repo $document) -Destination $stage
    }
    Compress-Archive -LiteralPath $stage -DestinationPath $zip -CompressionLevel Optimal
    $hash = (Get-FileHash -LiteralPath $zip -Algorithm SHA256).Hash.ToLowerInvariant()
    Set-Content -LiteralPath $sums -Value "$hash  $([System.IO.Path]::GetFileName($zip))" -Encoding ascii
    Write-Output "Archive: $zip"
    Write-Output "SHA256: $hash"
} finally {
    if (Test-Path -LiteralPath $stage) { Remove-Item -LiteralPath $stage -Recurse -Force }
    Pop-Location
}
