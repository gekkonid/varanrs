#!/usr/bin/env bash
# Fail if the release binary links any shared library other than the C
# runtime (we ship self-contained binaries).
set -euo pipefail
bin="${1:-target/release/varanrs}"
case "$(uname -s)" in
  Linux)
    deps=$(ldd "$bin" | awk '{print $1}')
    allowed='^(linux-vdso\.so\.1|linux-gate\.so\.1|libc\.so\.6|libm\.so\.6|libgcc_s\.so\.1|libpthread\.so\.0|libdl\.so\.2|librt\.so\.1|/lib.*/ld-linux.*\.so\.[0-9]+)$'
    ;;
  Darwin)
    deps=$(otool -L "$bin" | tail -n +2 | awk '{print $1}')
    allowed='^(/usr/lib/libSystem\.B\.dylib|/usr/lib/libiconv\.2\.dylib|/usr/lib/libc\+\+\.1\.dylib)$'
    ;;
  *) echo "unsupported OS"; exit 1 ;;
esac
echo "$deps"
bad=$(echo "$deps" | grep -Ev "$allowed" || true)
if [ -n "$bad" ]; then
  echo "ERROR: unexpected shared library dependencies:" >&2
  echo "$bad" >&2
  exit 1
fi
echo "linkage OK"
