#!/usr/bin/env bash
# Fetch Eigen 3.4.0 headers (required to compile shp_core).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EIGEN="$ROOT/third_party/eigen"
if [[ -f "$EIGEN/Eigen/Eigen" ]]; then
  echo "Eigen already present at $EIGEN"
  exit 0
fi
mkdir -p "$ROOT/third_party"
if command -v git >/dev/null 2>&1; then
  git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git "$EIGEN"
else
  echo "git is required to download Eigen." >&2
  exit 1
fi
echo "Eigen installed to $EIGEN"
