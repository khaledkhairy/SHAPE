#!/usr/bin/env bash
# Easy Ubuntu installer for SHAPE (Windows/Linux cross-platform port).
# Usage:
#   git clone https://github.com/khaledkhairy/SHAPE.git
#   cd SHAPE
#   ./install_ubuntu.sh
#
# After install, launch from the Applications menu ("SHAPE") or run:  shape
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

echo "==> Installing system packages (sudo may prompt)..."
apt-get update
apt-get install -y \
  python3 python3-venv python3-dev python3-pip \
  build-essential git \
  libgl1 libglib2.0-0 libxkbcommon-x11-0 libxcb-xinerama0 \
  libxcb-cursor0 libegl1 libxrender1 libfontconfig1 \
  libx11-xcb1 libdbus-1-3 libsm6 libice6

echo "==> Creating Python virtual environment..."
python3 -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

echo "==> Fetching Eigen headers..."
bash scripts/fetch_eigen.sh

echo "==> Building C++ core (shp_core)..."
python setup.py build_ext --inplace
SO="$(ls -1 shp_core*.so 2>/dev/null | head -1 || true)"
if [[ -z "$SO" ]]; then
  echo "ERROR: shp_core extension was not built." >&2
  exit 1
fi
echo "Built $SO"

LAUNCHER="$HOME/.local/bin/shape"
mkdir -p "$HOME/.local/bin"
cat > "$LAUNCHER" <<EOF
#!/usr/bin/env bash
cd "$ROOT/data" 2>/dev/null || cd "$ROOT"
exec "$ROOT/.venv/bin/python" "$ROOT/app/shape_app.py" "\$@"
EOF
chmod +x "$LAUNCHER"

DESKTOP="$HOME/.local/share/applications/shape.desktop"
mkdir -p "$(dirname "$DESKTOP")"
cat > "$DESKTOP" <<EOF
[Desktop Entry]
Type=Application
Name=SHAPE
Comment=Spherical Harmonics Parameterization Explorer
Exec=$LAUNCHER
Path=$ROOT
Terminal=false
Categories=Science;Education;
StartupWMClass=shape
EOF
chmod +x "$DESKTOP" 2>/dev/null || true

if [[ -d "$HOME/.local/share/icons/hicolor" ]]; then
  ICON="$HOME/.local/share/icons/hicolor/256x256/apps/shape.png"
  if [[ -f "$ROOT/clips/Screen Shot 2016-02-16 at 10.05.11 AM.png" ]]; then
    mkdir -p "$(dirname "$ICON")"
    cp "$ROOT/clips/Screen Shot 2016-02-16 at 10.05.11 AM.png" "$ICON" 2>/dev/null || true
    sed -i 's|^Exec=|Icon=shape\nExec=|' "$DESKTOP" 2>/dev/null || true
  fi
fi

echo ""
echo "=============================================="
echo " SHAPE installed successfully."
echo " Launch:  shape"
echo " Or find 'SHAPE' in your Applications menu."
echo "=============================================="
if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
  echo "Tip: add ~/.local/bin to your PATH if 'shape' is not found."
fi
