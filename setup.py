"""Build the shp_core extension module (pybind11 wrapper over the C++ core).

Usage:
    python setup.py build_ext --inplace

Eigen is vendored under third_party/eigen and the (nearly untouched) core
headers live under core/.
"""
import sys
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup


class BuildExtStaticCRT(build_ext):
    """On MSVC, statically link the C runtime (/MT) so the resulting .pyd has no
    external msvcp140.dll dependency. This avoids DLL-conflict crashes when the
    module is loaded alongside (or bundled with) PyQt5, which ships an older
    msvcp140.dll."""

    def build_extensions(self):
        if self.compiler.compiler_type == "msvc":
            for attr in ("compile_options", "compile_options_debug"):
                opts = getattr(self.compiler, attr, None)
                if not opts:
                    continue
                opts[:] = [o for o in opts if o.upper() not in ("/MD", "/MDD")]
            for ext in self.extensions:
                ext.extra_compile_args = (ext.extra_compile_args or []) + ["/MT"]
        super().build_extensions()

ROOT = Path(__file__).parent.resolve()
CORE = ROOT / "core"
EIGEN = ROOT / "third_party" / "eigen"

if sys.platform == "win32":
    extra_compile_args = ["/O2", "/EHsc", "/bigobj", "/D_CRT_SECURE_NO_WARNINGS", "/wd4996"]
    extra_link_args = []
else:
    extra_compile_args = ["-O3", "-w"]
    extra_link_args = []

ext_modules = [
    Pybind11Extension(
        "shp_core",
        [str(ROOT / "bindings" / "shp_core.cpp")],
        include_dirs=[str(CORE), str(EIGEN)],
        cxx_std=14,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

setup(
    name="shp_core",
    version="1.0.0",
    description="SPHARM surface core for the SHAPE viewer",
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExtStaticCRT},
)
