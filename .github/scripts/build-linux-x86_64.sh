#!/usr/bin/env bash
set -euo pipefail

export PATH="/root/.cargo/bin:${PATH}"

if ! command -v rustup >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain none
fi

rustup set profile minimal
rustup toolchain install "${RUST_TOOLCHAIN}"
rustup default "${RUST_TOOLCHAIN}"

yum install -y bcftools htslib htslib-tools

# bindgen (used by rust-htslib's hts-sys) needs libclang at build time only;
# take it from the PyPI wheel, which works on manylinux2014.
PY=/opt/python/cp312-cp312/bin/python
"$PY" -m pip install --quiet libclang
export LIBCLANG_PATH="$("$PY" -c 'import clang, os; print(os.path.join(os.path.dirname(clang.__file__), "native"))')"

cargo build --verbose --release
.github/scripts/check-linkage.sh target/release/varanrs
cargo test --verbose --release
