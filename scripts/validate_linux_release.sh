#!/usr/bin/env bash
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
version="$(sed -nE 's/^version[[:space:]]*=[[:space:]]*"([^"]+)"/\1/p' "$repo/Cargo.toml" | head -n 1)"
target="x86_64-unknown-linux-gnu"
dist="$repo/target/dist/linux"
archive="$dist/electrochem-tools-v${version}-${target}.tar.gz"
exe="$repo/target/$target/release/eiscli"

if [[ -z "$version" ]]; then
  echo "Could not read package version from Cargo.toml" >&2
  exit 1
fi

(
  cd "$dist"
  sha256sum --check SHA256SUMS-linux.txt
)

archive_entries="$(tar -tzf "$archive")"
for expected in eiscli clean_eis merge_cor trim_cv README.md LICENSE THIRD_PARTY_NOTICES.md CHANGELOG.md; do
  if ! grep -q "/${expected}$" <<< "$archive_entries"; then
    echo "Linux archive is missing $expected" >&2
    exit 1
  fi
done

smoke_root="$(mktemp -d /tmp/electrochem-tools-v013-smoke.XXXXXX)"

version_output="$("$exe" --version)"
if [[ "$version_output" != *"$version"* ]]; then
  echo "Unexpected version output: $version_output" >&2
  exit 1
fi

"$exe" --help >/dev/null
"$exe" validate --help >/dev/null
"$exe" clean --help >/dev/null
"$exe" drt --help >/dev/null
"$exe" fit-ecm --help >/dev/null
"$exe" validate -i "$repo/tests/fixtures/bayesian_eis.z60" \
  --n-tau 30 --kk-residual-threshold 3 --kk-min-points 12 \
  --out-root "$smoke_root/kk" >/dev/null
"$exe" clean -i "$repo/tests/fixtures/bayesian_eis.z60" \
  --out-root "$smoke_root/clean" >/dev/null
"$exe" drt -i "$repo/tests/fixtures/bayesian_eis.z60" \
  --nonnegative --n-tau 30 --kk-residual-threshold 3 --kk-min-points 12 \
  --out-root "$smoke_root/drt" >/dev/null
"$exe" fit-ecm -i "$repo/tests/fixtures/bayesian_eis.z60" \
  --model R_QR --auto-init --kk-residual-threshold 3 --kk-min-points 12 \
  --out-root "$smoke_root/ecm" >/dev/null

python3 - "$smoke_root" <<'PY'
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
kk = json.loads((root / "kk/bayesian_eis_kk_summary.json").read_text())
drt_trim = json.loads(
    (root / "drt/bayesian_eis_drt/kk_trim_summary.json").read_text()
)
drt_summary = json.loads(
    (root / "drt/bayesian_eis_drt/residual_summary.json").read_text()
)
ecm_trim = json.loads(
    (root / "ecm/bayesian_eis_fit_ecm/kk_trim_summary.json").read_text()
)
ecm_summary = json.loads(
    (root / "ecm/bayesian_eis_fit_ecm/fit_params.json").read_text()
)

assert kk["trimmed_high_frequency_points"] >= 1
assert drt_trim["retained_points"] == drt_summary["n_points"]
assert ecm_trim["retained_points"] == ecm_summary["n_points"]
PY

echo "Linux release assets, checksums, archive, and smoke tests: OK ($version_output)"
