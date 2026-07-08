# SHAPE (Spherical HArmonics Parameterization Explorer)

![Alt text](https://github.com/khaledkhairy/SHAPE/blob/master/clips/Screen%20Shot%202016-02-16%20at%2010.05.11%20AM.png "SHAPE screenshot")
-----------------------------------------------------------------------------
Copyright (c) 2017, HHMI-Janelia Research Campus All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  
* Neither the name HHMI-Janelia Research Campus nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL HHMI-JANELIA RESEARCH CAMPUS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Status: 
------------------------------------------------------------------------------
In production use at Janelia. This is a nascent set of tools that is undergoing large changes and code cleanup. We consider the library suitable for use by our collaborators as well as other research groups. Due to limited staffing, we do not guarantee support for outside groups.
## Summary: 
------------------------------------------------------------------------------
Generate spherical harmonics-based shapes, or explore existing ones, and their shape properties interactively and accurately.
SHAPE allows monitoring behavior of (and modification of) spherical harmonics parameterization (SPHARM) coefficients interactively.
These are (Fourier) coefficients of a parametric spherical harmonics representation of objects of genus zero.
Genus zero objects have no handles (the topology is that of the sphere).
SHAPE is also used to view spherical harmonics fields mapped onto a base SPHARM shape.
Shapes generated with SHAPE can be exported as triangular surface meshes into OBJ, PLY, STL and VRML formats.
Test shapes are provided for initial experimentation.

For examples of published works based on shape parameterization using Fourier series expansions, please refer to:

-Minimum-energy vesicle and cell shapes calculated using spherical harmonics parameterization. 
Khaled Khairy and Jonathon Howard. Soft Matter, 2011,7, 2138-2143

-Drawing an elephant with 4 parameters. 
Mayer J., K. Khairy and J. Howard, American Journal of Physics, 78(6):648-649 (2010).

-Shapes of Red Blood Cells: Comparison of 3D Confocal Images with the Bilayer-Couple Model. 
Khaled Khairy, JiJinn Foo and Jonathon Howard. Cellular and Molecular Bioengineering, 1:173-181, 2008

-Detection of Deformable Objects in 3D Images using Markov Chain Monte Carlo and Spherical Harmonics. 
Khairy K., E. Reynaud, E. Stelzer, MICCAI New York 2008

-Spherical harmonics-based parametric deconvolution of 3D surface images using bending energy minimization. 
Khairy Khairy and Jonathon Howard. Med Image Anal. 2008 Apr;12(2):217-27. Epub 2007 Oct 30.

-Spherical Harmonics 3D Active Contours for Membrane Bilayer-Bound Surfaces. 
Khaled Khairy, Jacques Pecreaux and Jonathon Howard. MICCAI 2007



Building notes:
------------------------------------------------------------------------------

SHAPE has been built using the following C/C++ libraries on both MacOSX and Windows 7:

[1] Eigen 3.1.1: http://eigen.tuxfamily.org

[2] VTK 7.1.1: http://www.vtk.org/

[3] Qt 5.7: http://www.qt.io/

[4] shape tools: https://github.com/khaledkhairy/shape_tools

Building has not been tested extensively on any other configuration.
For MacOSX Sierra, it is recommended to use the supplied binary.

SHAPE was compiled using Qt 5.7and VTK 7.1.1 (although other versions might also work)

## Acknowledgements:
------------------------------------------------------------------------------
To Jonathon Howard for supporting development of the first version of SHAPE.

---

## Windows / Linux port (2026)

The original macOS C++ application (`main.cpp`, `shape.cpp`, `SHAPE_binary_MACOSX/`, etc.) is **unchanged**. A cross-platform port reuses the same `shape_tools` math core and the same GUI layout, hosted through PyQt5 + VTK:

| Component | Location |
|---|---|
| Original macOS app | repo root (`main.cpp`, `shape.cpp`, `CMakeLists.txt`, `SHAPE_binary_MACOSX/`) |
| Cross-platform core | `core/` (vendored `shape_tools`; only Eigen include path repointed) |
| Cross-platform GUI | `app/shape.ui`, `app/shape_app.py` |
| pybind11 wrapper | `bindings/shp_core.cpp`, `setup.py` |

Eigen headers are downloaded automatically by `scripts/fetch_eigen.sh` (not stored in git).

---

## Installation — Windows

Tested on **Windows 10/11 (64-bit)**.

### Prerequisites

1. **Python 3.10 or newer (64-bit)**  
   Download from [python.org](https://www.python.org/downloads/). During setup, check **“Add python.exe to PATH”**.

2. **Microsoft C++ Build Tools (MSVC)**  
   Required to compile the `shp_core` extension (the C++ math core). Install via:

   ```powershell
   winget install Microsoft.VisualStudio.2022.BuildTools
   ```

   When prompted, include the **“Desktop development with C++”** workload (MSVC v143 + Windows 10/11 SDK). Accept any UAC elevation prompt.

3. **Git** (optional but recommended)  
   To clone the repository:

   ```powershell
   winget install Git.Git
   ```

### Option A — Build the double-clickable application (recommended)

Open **PowerShell**, clone the repo, and run the build script:

```powershell
git clone https://github.com/khaledkhairy/SHAPE.git
cd SHAPE
.\build.ps1
```

`build.ps1` performs these steps automatically:

1. Creates a Python virtual environment (`.venv/`)
2. Installs dependencies from `requirements.txt` (PyQt5, VTK, pybind11, PyInstaller, …)
3. Downloads Eigen headers and compiles `shp_core` with MSVC
4. Packages a standalone folder with PyInstaller

When finished, launch the app by double-clicking:

```
dist\SHAPE\SHAPE.exe
```

You can copy the entire `dist\SHAPE\` folder anywhere on your machine; it is self-contained.

**Open a shape file from the command line:**

```powershell
.\dist\SHAPE\SHAPE.exe "C:\path\to\your_shape.shp3"
```

Sample shapes are included under `data\test_shapes\`.

### Option B — Run from source (for development)

Faster iteration when editing code; does not produce a packaged `.exe`:

```powershell
git clone https://github.com/khaledkhairy/SHAPE.git
cd SHAPE
python -m venv .venv
.\.venv\Scripts\pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pip install pybind11
bash scripts/fetch_eigen.sh          # requires Git Bash, or fetch Eigen manually (see below)
.\.venv\Scripts\python.exe setup.py build_ext --inplace
.\run_from_source.ps1
```

Or, after the first build, open a specific shape:

```powershell
.\run_from_source.ps1 data\test_shapes\RBC_and_Vesicles\discocyte.shp3
```

### Manual build steps (Windows)

If you prefer not to use the helper scripts:

```powershell
# 1. Virtual environment
python -m venv .venv
.\.venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt

# 2. Eigen headers (requires git)
git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git third_party/eigen

# 3. Compile C++ core
python setup.py build_ext --inplace
# Produces: shp_core.cp312-win_amd64.pyd  (version tag varies with Python)

# 4a. Run directly
python app\shape_app.py

# 4b. Package as double-clickable app
pyinstaller shape.spec --noconfirm --distpath dist --workpath build_pyi
# Output: dist\SHAPE\SHAPE.exe
```

### Windows troubleshooting

| Problem | Solution |
|---|---|
| `cl` / compiler not found when building | Open a **“x64 Native Tools Command Prompt for VS 2022”** or reinstall Build Tools with the C++ workload. `build.ps1` normally finds MSVC automatically via setuptools. |
| App crashes immediately on launch | Rebuild `shp_core` with `python setup.py build_ext --inplace` after any Python or MSVC update. The extension is built with a static C runtime (`/MT`) to avoid conflicts with PyQt5’s bundled DLLs. |
| PyInstaller build fails (locked files) | Close any running `SHAPE.exe` instances, then rerun `.\build.ps1`. |
| Drag-and-drop / Open does nothing | Ensure the file has a `.shp3` extension and is a valid SHAPE coefficient file. |

---

## Installation — Linux (Ubuntu)

Tested on **Ubuntu 22.04 LTS** and **Ubuntu 24.04 LTS**. Other Debian-based distributions should work with the same package names.

> **Important:** clone the repository onto a **native Linux filesystem** (ext4, btrfs, etc.). Do not build inside a Windows-mounted drive (Dropbox, `/mnt/c/`, …) — file permissions on those mounts can break the Python virtual environment.

### Option A — One-step install (recommended)

```bash
git clone https://github.com/khaledkhairy/SHAPE.git
cd SHAPE
chmod +x install_ubuntu.sh
./install_ubuntu.sh
```

The installer will:

1. Install system packages (`python3`, `build-essential`, OpenGL/X11 libraries) via `apt`
2. Create `.venv/` and install Python dependencies
3. Download Eigen 3.4.0 headers into `third_party/eigen/`
4. Compile `shp_core.so` (the C++ math core)
5. Install a `shape` command in `~/.local/bin/`
6. Register **SHAPE** in your desktop Applications menu

Launch the app:

```bash
shape
```

Or find **SHAPE** in the Applications menu (Science / Education category).

Open a shape file directly:

```bash
shape ~/SHAPE/data/test_shapes/RBC_and_Vesicles/discocyte.shp3
```

If `shape` is not found, add `~/.local/bin` to your PATH:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Option B — Run from source (for development)

After a one-time install (`./install_ubuntu.sh`), rebuild and launch with:

```bash
./build_linux.sh
```

Or pass a shape file:

```bash
./build_linux.sh data/test_shapes/General_created/sphere_default.shp3
```

### Manual build steps (Linux)

```bash
# 1. System dependencies
sudo apt-get update
sudo apt-get install -y \
  python3 python3-venv python3-dev python3-pip \
  build-essential git \
  libgl1 libglib2.0-0 libxkbcommon-x11-0 libxcb-xinerama0 \
  libxcb-cursor0 libegl1 libxrender1 libfontconfig1 \
  libx11-xcb1 libdbus-1-3 libsm6 libice6

# 2. Clone and enter repo
git clone https://github.com/khaledkhairy/SHAPE.git
cd SHAPE

# 3. Python environment
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

# 4. Eigen headers
bash scripts/fetch_eigen.sh

# 5. Compile C++ core
python setup.py build_ext --inplace
# Produces: shp_core.cpython-312-x86_64-linux-gnu.so  (tag varies with Python)

# 6. Run
cd data    # mesh .tri tables are read from the working directory
python ../app/shape_app.py
```

### Linux troubleshooting

| Problem | Solution |
|---|---|
| `shape: command not found` | Ensure `~/.local/bin` is on your PATH (see above), or run `~/.local/bin/shape`. |
| `ImportError: libGL.so.1` or X11 errors | Reinstall graphics libraries: `sudo apt-get install -y libgl1 libxcb-cursor0 libegl1`. |
| `shp_core` compile error: Eigen not found | Run `bash scripts/fetch_eigen.sh` and confirm `third_party/eigen/Eigen/Eigen` exists. |
| Permission errors during `pip install` | Do not run the installer as root. Clone to a directory you own on a native Linux filesystem. |
| Mesh-approx curvature unavailable | Copy the appropriate `.tri` file (`uni900.tri`, `uni10k.tri`, …) into the working directory, or launch via `shape` / `install_ubuntu.sh` which sets the working directory to `data/`. |

### Other Linux distributions

The Python/VTK/PyQt5 stack is distribution-agnostic once system libraries are present. On Fedora/RHEL:

```bash
sudo dnf install python3 python3-devel gcc-c++ git mesa-libGL
# then follow the manual build steps above
```

---

## Using SHAPE (Windows and Linux)

After installation, the app opens on a **unit sphere**. From there:

- **Sliders** (49 per axis × X, Y, Z) edit SPHARM coefficients; the surface updates in real time.
- **Drag and drop** a `.shp3` file onto the window, or use **File → Open**.
- **File → Save** writes coefficients back to `.shp3`.
- **File → Export** saves the current mesh as OBJ, STL, PLY, or VRML.
- The **scalar-field** combo box colours the surface (θ, φ, mean/Gaussian curvature, coordinates, or fields stored in the file).
- **Lmax**, mesh resolution, and Gaussian-quadrature density can be changed from the control panel; area, volume, reduced volume, and bending energy update live.

Sample shapes are in `data/test_shapes/` (red blood cells, vesicles, demo shapes, etc.).
