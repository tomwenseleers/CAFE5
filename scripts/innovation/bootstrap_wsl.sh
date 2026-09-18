#!/usr/bin/env bash
# Optional user-local tooling when CMake/SciPy are unavailable and sudo is not used.
# Ubuntu package names: verified on Ubuntu resolute under WSL2.
set -euo pipefail
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
tools_dir="$repo_dir/.wsl-tools"
mkdir -p "$tools_dir/packages" "$tools_dir/root"
cd "$tools_dir/packages"
apt-get download cmake cmake-data librhash1 libjsoncpp26 libuv1t64 \
  libarchive13t64 libcurl4t64 libexpat1 libqhull-r8.0 \
  python3-numpy python3-scipy python3-decorator
for package in ./*.deb; do dpkg-deb -x "$package" "$tools_dir/root"; done
cat > "$tools_dir/env.sh" <<'ENV'
cafe_tools_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export PATH="$cafe_tools_dir/root/usr/bin:$PATH"
export LD_LIBRARY_PATH="$cafe_tools_dir/root/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="$cafe_tools_dir/root/usr/lib/python3/dist-packages${PYTHONPATH:+:$PYTHONPATH}"
unset cafe_tools_dir
ENV
printf 'Tooling ready. Run: source "%s/env.sh"\n' "$tools_dir"
