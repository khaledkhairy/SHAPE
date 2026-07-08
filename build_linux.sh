#!/usr/bin/env bash
# Rebuild shp_core and run SHAPE from source on Linux.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
if [[ ! -d .venv ]]; then
  echo "Run ./install_ubuntu.sh first." >&2
  exit 1
fi
# shellcheck disable=SC1091
source .venv/bin/activate
bash scripts/fetch_eigen.sh
python setup.py build_ext --inplace
exec python app/shape_app.py "$@"
