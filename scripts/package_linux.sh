#!/usr/bin/env bash
set -euo pipefail

version="${1:-}"
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
manifest="$repo/Cargo.toml"
if [[ -z "$version" ]]; then
  version="$(sed -nE 's/^version[[:space:]]*=[[:space:]]*"([^"]+)"/\1/p' "$manifest" | head -n 1)"
fi
if [[ -z "$version" ]]; then
  echo "Could not read package version from Cargo.toml" >&2
  exit 1
fi

target="x86_64-unknown-linux-gnu"
name="electrochem-tools-v${version}-${target}"
dist="$repo/target/dist/linux"
stage="$repo/target/package-stage/$name"
archive="$dist/$name.tar.gz"
sums="$dist/SHA256SUMS-linux.txt"
binaries=(eiscli clean_eis merge_cor trim_cv)
cargo_config=()

if [[ -n "${CARGO_VENDOR_DIR:-}" ]]; then
  if [[ ! -d "$CARGO_VENDOR_DIR" ]]; then
    echo "Cargo vendor directory does not exist: $CARGO_VENDOR_DIR" >&2
    exit 1
  fi
  cargo_config+=(
    --config 'source.crates-io.replace-with="vendored-sources"'
    --config "source.vendored-sources.directory=\"$CARGO_VENDOR_DIR\""
  )
fi

case "$dist" in "$repo"/*) ;; *) echo "Refusing dist path outside repository: $dist" >&2; exit 1 ;; esac
case "$stage" in "$repo"/*) ;; *) echo "Refusing stage path outside repository: $stage" >&2; exit 1 ;; esac

rm -rf -- "$dist" "$stage"
mkdir -p -- "$dist" "$stage"

cargo "${cargo_config[@]}" build --locked --release --all-features --target "$target" --jobs 1

assets=()
for binary in "${binaries[@]}"; do
  source="$repo/target/$target/release/$binary"
  if [[ ! -f "$source" ]]; then
    echo "Missing binary: $source" >&2
    exit 1
  fi
  install -m 755 "$source" "$stage/$binary"
  standalone="$dist/${binary}-linux-x86_64"
  install -m 755 "$source" "$standalone"
  assets+=("$standalone")
done

for document in README.md LICENSE THIRD_PARTY_NOTICES.md CHANGELOG.md; do
  install -m 644 "$repo/$document" "$stage/$document"
done

tar -C "$(dirname "$stage")" -czf "$archive" "$name"
assets+=("$archive")

: > "$sums"
for asset in "${assets[@]}"; do
  hash="$(sha256sum "$asset" | awk '{print $1}')"
  printf '%s  %s\n' "$hash" "$(basename "$asset")" >> "$sums"
done

printf 'Archive: %s\n' "$archive"
printf 'Standalone binaries: %s\n' "${binaries[*]}"
printf 'Checksums: %s\n' "$sums"

rm -rf -- "$stage"
