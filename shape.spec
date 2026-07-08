# PyInstaller spec for SHAPE (Spherical HArmonics Parameterization Explorer)
# Build with:  .\.venv\Scripts\pyinstaller.exe shape.spec --noconfirm
import glob
import os

from PyInstaller.utils.hooks import collect_submodules

block_cipher = None

# The compiled C++ core (built by setup.py build_ext --inplace)
pyd = glob.glob("shp_core*.pyd")
assert pyd, "shp_core*.pyd not found - run: python setup.py build_ext --inplace"

binaries = [(pyd[0], ".")]

datas = [
    ("app/shape.ui", "."),
    ("data", "data"),
]

hiddenimports = [
    "shp_core",
    "vtkmodules.qt.QVTKRenderWindowInteractor",
    "vtkmodules.util.numpy_support",
    "vtkmodules.all",
] + collect_submodules("vtkmodules")

a = Analysis(
    ["app/shape_app.py"],
    pathex=["."],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=["matplotlib", "tkinter", "pyvista", "numba"],
    cipher=block_cipher,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="SHAPE",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="SHAPE",
)
