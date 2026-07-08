"""SHAPE - Spherical HArmonics Parameterization Explorer (Windows/Linux port).

This is a faithful port of the original SHAPE viewer glue logic (shape.cpp) to
PyQt5 + VTK. It loads the *original* Qt Designer layout (shape.ui) and drives the
*original* C++ core math (compiled as the `shp_core` pybind11 module) so the
numerical behaviour is identical to the original application.

Copyright (c) 2017, HHMI-Janelia Research Campus (original SHAPE).
"""
import os
import sys

# ---------------------------------------------------------------------------
# Resource resolution (works both in-source and when frozen by PyInstaller)
# ---------------------------------------------------------------------------
def resource_dir() -> str:
    if getattr(sys, "frozen", False):
        base = getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
    else:
        base = os.path.dirname(os.path.abspath(__file__))
    return base


APP_DIR = resource_dir()
# Make the compiled core importable when running from source.
sys.path.insert(0, os.path.dirname(APP_DIR))
sys.path.insert(0, APP_DIR)

# The C++ core reads spherical-mesh connectivity files (e.g. uni900.tri) from the
# current working directory. Point CWD at the bundled data folder so mesh-based
# curvature works out of the box.
DATA_DIR = os.path.join(APP_DIR, "data")
if not os.path.isdir(DATA_DIR):
    # running from source layout: <root>/data
    DATA_DIR = os.path.join(os.path.dirname(APP_DIR), "data")
if os.path.isdir(DATA_DIR):
    try:
        os.chdir(DATA_DIR)
    except OSError:
        pass

# IMPORTANT: import the compiled core BEFORE PyQt5/VTK. shp_core is built with
# the VS2022 C runtime; PyQt5 ships an older msvcp140.dll and, if loaded first,
# would satisfy shp_core's runtime with an incompatible DLL and crash on load.
# Importing shp_core first pins the newer (backward-compatible) runtime.
import shp_core

import numpy as np
from PyQt5 import QtCore, QtWidgets, uic
from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow

import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.util import numpy_support

UI_PATH = os.path.join(APP_DIR, "shape.ui")


class ShapeWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        uic.loadUi(UI_PATH, self)
        self.setAcceptDrops(True)

        # ---- state (mirrors shape.cpp members) ----------------------------
        self._minH = 0.0
        self._maxH = 0.0
        self._scale_fac = 3.0
        self._axis_on = True
        self.scalar_bar_vis = True
        self.cube_axis_actor1_on = False
        self.check_self_intersection = True
        self.cscale = 1.0

        # ---- VTK render window embedded into the vtkWidget1 placeholder ----
        self.vtk_iren = QVTKRenderWindowInteractor(self.vtkWidget1)
        layout = QtWidgets.QVBoxLayout(self.vtkWidget1)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.vtk_iren)
        self._ren1 = vtk.vtkRenderer()
        self.vtk_iren.GetRenderWindow().AddRenderer(self._ren1)
        self._ren1.SetBackground(0.0, 0.0, 0.0)

        # ---- core surface object (Lmax=24, gauss dim=30, mesh tri_n=4) -----
        self.shp = shp_core.Shape(24, 30, 4)
        self.cscale = 1.0 / self._cscale_calc()

        # ---- VTK pipeline objects -----------------------------------------
        self.lut = vtk.vtkLookupTable()
        self.lut.SetNumberOfTableValues(400)
        self.lut.SetHueRange(0.667, 0.0)
        self.lut.SetSaturationRange(1, 1)
        self.lut.SetValueRange(1, 1)
        self.lut.SetRampToSCurve()

        self.polydata = vtk.vtkPolyData()
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetLookupTable(self.lut)
        self.mapper.SetInputData(self.polydata)
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self._ren1.AddActor(self.actor)

        # axes actor
        self._axes_actor = vtk.vtkAxesActor()
        self._axes_actor.SetShaftTypeToLine()
        self._axes_actor.SetConeRadius(0.1)
        self._axes_actor.SetXAxisLabelText("x")
        self._axes_actor.SetYAxisLabelText("y")
        self._axes_actor.SetZAxisLabelText("z")
        self._axes_actor.SetTotalLength(5, 5, 5)
        self._ren1.AddActor(self._axes_actor)

        # cube axes actor (toggled)
        self._cube_axes_actor1 = vtk.vtkCubeAxesActor()

        # scalar bar
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(1, 1, 1)
        tprop.SetFontFamilyToArial()
        tprop.SetFontSize(32)
        self.scalarBar = vtk.vtkScalarBarActor()
        self.scalarBar.SetNumberOfLabels(4)
        self.scalarBar.SetMaximumWidthInPixels(80)
        self.scalarBar.SetLabelTextProperty(tprop)
        self.scalarBar.SetTitleTextProperty(tprop)
        self.scalarBar.SetVisibility(False)
        self._ren1.AddActor(self.scalarBar)

        # face connectivity cache
        self._cells = None
        self._rebuild_cells()

        # ---- wire up widgets ----------------------------------------------
        self._connect_signals()

        # ---- initial state (mirrors constructor tail) ---------------------
        self.on_initializeButton_clicked()

        self.comboBox_gdim.blockSignals(True)
        self.comboBox_gdim.setCurrentIndex(4)          # dim = 30
        self.comboBox_gdim.blockSignals(False)
        self.comboBox_spherical_mesh.blockSignals(True)
        self.comboBox_spherical_mesh.setCurrentIndex(4 - 1)  # tri_n = 4
        self.comboBox_spherical_mesh.blockSignals(False)
        self.checkBox_self_intersection.blockSignals(True)
        self.checkBox_self_intersection.setChecked(True)
        self.checkBox_self_intersection.blockSignals(False)
        # self-intersection is O(n^2): disabled above tri_n 2 (as in original)
        self.checkBox_self_intersection.setChecked(False)
        self.checkBox_self_intersection.setEnabled(False)
        self.check_self_intersection = False

        self.setWindowTitle("SHAPE v. 1.0")
        self.vtk_iren.Initialize()
        self._ren1.ResetCamera()
        self.vtk_iren.GetRenderWindow().Render()

    # ======================================================================
    # signal wiring
    # ======================================================================
    def _connect_signals(self):
        self.initializeButton.clicked.connect(self.on_initializeButton_clicked)
        self.axes_toggle_Button.clicked.connect(self.on_axes_toggle_Button_clicked)
        self.pushButton_cube_axis_toggle.clicked.connect(
            self.on_pushButton_cube_axis_toggle_clicked)
        self.lineEdit_Lmax.returnPressed.connect(self.on_lineEdit_Lmax_returnPressed)

        self.comboBox_gdim.currentIndexChanged.connect(self.on_comboBox_gdim_changed)
        self.comboBox_spherical_mesh.currentIndexChanged.connect(
            self.on_comboBox_spherical_mesh_changed)
        self.comboBox_sf.currentIndexChanged.connect(self.on_comboBox_sf_changed)
        self.checkBox_self_intersection.stateChanged.connect(
            self.on_checkBox_self_intersection_changed)

        self.actionOpen.triggered.connect(self.on_actionOpen)
        self.actionSave.triggered.connect(self.on_actionSave)
        self.actionExport.triggered.connect(self.on_actionExport)
        self.actionExport_to_STL.triggered.connect(self.on_actionExport_to_STL)
        self.actionExport_to_PLY.triggered.connect(self.on_actionExport_to_PLY)
        self.actionExport_to_VRML.triggered.connect(self.on_actionExport_to_VRML)
        self.actionExit.triggered.connect(self.on_actionExit)

        # 49 coefficients x 3 axes = 147 sliders
        # coordix: x=0, y=1, z=2 (as in shape.cpp)
        self._sliders = []
        for i in range(1, 50):
            cix = i - 1
            for suffix, coordix, setter in (
                ("", 2, self.shp.set_zc),      # z
                ("xx", 0, self.shp.set_xc),    # x
                ("yy", 1, self.shp.set_yc),    # y
            ):
                name = "verticalSliderx%02d%s" % (i, suffix)
                w = getattr(self, name, None)
                if w is None:
                    continue
                self._sliders.append((w, coordix, cix, setter))
                w.valueChanged.connect(
                    lambda val, cx=coordix, ci=cix, st=setter: self._on_slider(cx, ci, st, val))

    # ======================================================================
    # slider handling (mirrors on_verticalSliderXX_valueChanged)
    # ======================================================================
    def _get_coeff(self, coordix, cix):
        if coordix == 0:
            return self.shp.get_xc(cix)
        if coordix == 1:
            return self.shp.get_yc(cix)
        return self.shp.get_zc(cix)

    def _on_slider(self, coordix, cix, setter, val):
        old_val = self._get_coeff(coordix, cix)
        new_val = float(val - 49) / 49.0 / self.cscale
        setter(cix, new_val)
        self.shp.set_needs_updating(True)
        self._update_vtkwindow_incremental(coordix, cix, old_val, new_val)

    # ======================================================================
    # core <-> vtk synchronisation
    # ======================================================================
    def _cscale_calc(self):
        mc = 0.0
        n = self.shp.n_coeffs()
        for i in range(1, min(4, n)):
            mc = max(mc, abs(self.shp.get_xc(i)), abs(self.shp.get_yc(i)),
                     abs(self.shp.get_zc(i)))
        return mc if mc > 1e-12 else 1.0

    def _rebuild_cells(self):
        faces = self.shp.get_faces().astype(np.int64)
        m = faces.shape[0]
        conn = np.empty((m, 4), dtype=np.int64)
        conn[:, 0] = 3
        conn[:, 1:] = faces
        idarr = numpy_support.numpy_to_vtkIdTypeArray(conn.ravel(), deep=1)
        cells = vtk.vtkCellArray()
        cells.SetCells(m, idarr)
        self._cells = cells

    def _update_vtkwindow_full(self):
        self.shp.update_full()
        self._update_gui()

    def _update_vtkwindow_incremental(self, coordix, cix, old_val, new_val):
        self.shp.update_incremental(coordix, cix, old_val, new_val)
        self._update_gui()

    def _build_polydata(self):
        pts = self.shp.get_points()
        vtk_pts = vtk.vtkPoints()
        vtk_pts.SetData(numpy_support.numpy_to_vtk(np.ascontiguousarray(pts), deep=1))
        self.polydata.SetPoints(vtk_pts)
        self.polydata.SetPolys(self._cells)

        sf_index = self.comboBox_sf.currentIndex()
        self._minH = 0.0
        self._maxH = 0.0
        if sf_index:
            sf = self.shp.compute_scalar_field(sf_index)
            self._minH = min(0.0, float(sf.min()))
            self._maxH = max(0.0, float(sf.max()))
            scalars = numpy_support.numpy_to_vtk(np.ascontiguousarray(sf), deep=1)
            self.polydata.GetPointData().SetScalars(scalars)
        else:
            self.polydata.GetPointData().SetScalars(None)

    def _update_gui(self):
        self._build_polydata()

        if self.check_self_intersection:
            if self.shp.self_intersect():
                print("Self intersection detected")

        # readouts
        self.lineEdit_Lmax.setText(str(self.shp.L_max()))
        self.lcdNumberArea.display(self.shp.A())
        self.lcdNumberVolume.display(self.shp.Vol())
        self.lcdNumberReducedVolume.display(self.shp.reduced_v())
        self.lcdNumberEb.display(self.shp.Eb())
        self.lcdNumberh.display(self.shp.h())
        self.lcdNumberWb.display(self.shp.T())
        self.lcdNumberArea_2.display(self.shp.smA())
        self.lcdNumberVolume_2.display(self.shp.smV())
        self.lcdNumberReducedVolume_2.display(self.shp.smv())
        self.lcdNumberEb_2.display(self.shp.smEb())
        self.lcdNumberh_2.display(self.shp.smh())

        sf_index = self.comboBox_sf.currentIndex()
        self.mapper.SetInputData(self.polydata)
        self.mapper.SetScalarRange(self._minH, self._maxH)
        if sf_index:
            self.mapper.ScalarVisibilityOn()
            self.mapper.SetScalarModeToUsePointData()
            self.mapper.SetColorModeToMapScalars()
        else:
            self.mapper.ScalarVisibilityOff()

        if self.cube_axis_actor1_on:
            self._cube_axes_actor1.SetBounds(self.actor.GetBounds())
            self._cube_axes_actor1.SetCamera(self._ren1.GetActiveCamera())
            self._ren1.AddActor(self._cube_axes_actor1)
        else:
            self._ren1.RemoveActor(self._cube_axes_actor1)

        if self.scalar_bar_vis and sf_index:
            scalars = self.polydata.GetPointData().GetScalars()
            if scalars is not None:
                rng = scalars.GetRange()
                self.lut.SetRange(rng)
                self.lut.Build()
                self.scalarBar.SetLookupTable(self.mapper.GetLookupTable())
                self.scalarBar.SetVisibility(True)
        else:
            self.scalarBar.SetVisibility(False)

        self.vtk_iren.GetRenderWindow().Render()

    def _synchronize_shape2sliders(self):
        for w, coordix, cix, _setter in self._sliders:
            val = self._get_coeff(coordix, cix)
            w.blockSignals(True)
            w.setValue(int(val * 49 * self.cscale) + 49)
            w.blockSignals(False)

    def _update_comboBox_sf(self):
        tags = self.shp.sf_tags()
        if not tags:
            return
        self.comboBox_sf.blockSignals(True)
        self.comboBox_sf.clear()
        for item in ("none", "theta", "phi",
                     "mean curvature (analytical)",
                     "mean curvature (mesh approx.)",
                     "Gaussian curvature (analytical)"):
            self.comboBox_sf.addItem(item)
        for t in tags:
            self.comboBox_sf.addItem(t)
        self.comboBox_sf.blockSignals(False)

    # ======================================================================
    # buttons / actions
    # ======================================================================
    def on_initializeButton_clicked(self):
        self.shp.initialize_to_sphere()
        self.cscale = 1.0 / self._cscale_calc()
        self._axes_actor.SetTotalLength(self._scale_fac / self.cscale,
                                        self._scale_fac / self.cscale,
                                        self._scale_fac / self.cscale)
        self.comboBox_sf.blockSignals(True)
        self.comboBox_sf.setCurrentIndex(0)
        self.comboBox_sf.blockSignals(False)
        self._synchronize_shape2sliders()
        self._update_vtkwindow_full()
        self._update_comboBox_sf()
        self._ren1.ResetCamera()
        self.vtk_iren.GetRenderWindow().Render()

    def on_axes_toggle_Button_clicked(self):
        self._axis_on = not self._axis_on
        self._axes_actor.SetVisibility(self._axis_on)
        self.vtk_iren.GetRenderWindow().Render()

    def on_pushButton_cube_axis_toggle_clicked(self):
        self.cube_axis_actor1_on = not self.cube_axis_actor1_on
        self._update_vtkwindow_full()

    def on_lineEdit_Lmax_returnPressed(self):
        try:
            Lmax_new = int(self.lineEdit_Lmax.text())
        except ValueError:
            Lmax_new = self.shp.L_max()
        Lmax_new = max(6, min(80, Lmax_new))
        self.shp.set_L_max(Lmax_new)
        self.cscale = 1.0 / self._cscale_calc()
        self._axes_actor.SetTotalLength(self._scale_fac / self.cscale,
                                        self._scale_fac / self.cscale,
                                        self._scale_fac / self.cscale)
        self._synchronize_shape2sliders()
        self._update_vtkwindow_full()

    def on_comboBox_gdim_changed(self, val):
        table = {0: 5, 1: 10, 2: 15, 3: 20, 4: 30, 5: 40, 6: 50,
                 7: 60, 8: 80, 9: 120, 10: 240, 11: 320}
        new_dim = table.get(val, 20)
        self.shp.set_new_basis(self.shp.L_max(), new_dim)
        self.shp.set_needs_updating(True)
        self._update_vtkwindow_full()

    def on_comboBox_spherical_mesh_changed(self, val):
        new_n = val + 1  # 0..5 -> 1..6
        new_n = max(1, min(6, new_n))
        self.shp.set_new_spherical_mesh(self.shp.L_max(), new_n)
        self._rebuild_cells()
        if new_n > 2:
            self.checkBox_self_intersection.setChecked(False)
            self.checkBox_self_intersection.setEnabled(False)
            self.check_self_intersection = False
        else:
            self.checkBox_self_intersection.setEnabled(True)
            if self.checkBox_self_intersection.isChecked():
                self.check_self_intersection = True
        self.shp.set_needs_updating(True)
        self._update_vtkwindow_full()

    def on_comboBox_sf_changed(self, val):
        self._update_vtkwindow_full()

    def on_checkBox_self_intersection_changed(self, val):
        self.check_self_intersection = bool(val)
        self._update_vtkwindow_full()

    # ---- file IO ----------------------------------------------------------
    def _shapes_dir(self):
        d = os.path.join(DATA_DIR, "test_shapes")
        return d if os.path.isdir(d) else DATA_DIR

    def on_actionOpen(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open shape", self._shapes_dir(), "SHP Files (*.shp3)")
        if filename:
            print("Opening:", filename)
            self.import_shape_from_disc(filename)

    def on_actionSave(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save shape to file ...", os.path.expanduser("~"),
            "SHP Files (*.shp3)")
        if filename:
            if not filename.lower().endswith(".shp3"):
                filename += ".shp3"
            self.shp.write(filename)

    def on_actionExport(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export to OBJ", os.path.expanduser("~"), "OBJ (*.obj)")
        if not filename:
            return
        prefix = filename[:-4] if filename.lower().endswith(".obj") else filename
        rw = self.vtk_iren.GetRenderWindow()
        axes_was_on = self._axis_on
        if axes_was_on:
            self._axes_actor.SetVisibility(False)
        rw.Render()
        obj = vtk.vtkOBJExporter()
        obj.SetRenderWindow(rw)
        obj.SetFilePrefix(prefix)
        obj.Write()
        if axes_was_on:
            self._axes_actor.SetVisibility(True)
        rw.Render()

    def on_actionExport_to_STL(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export to STL", os.path.expanduser("~"), "STL (*.stl)")
        if not filename:
            return
        if not filename.lower().endswith(".stl"):
            filename += ".stl"
        w = vtk.vtkSTLWriter()
        w.SetInputData(self.polydata)
        w.SetFileName(filename)
        w.Write()

    def on_actionExport_to_PLY(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export to PLY", os.path.expanduser("~"), "PLY (*.ply)")
        if not filename:
            return
        if not filename.lower().endswith(".ply"):
            filename += ".ply"
        ply = vtk.vtkPLYWriter()
        ply.SetInputData(self.polydata)
        ply.SetFileName(filename)
        ply.Write()

    def on_actionExport_to_VRML(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export to VRML", os.path.expanduser("~"), "VRML (*.wrl)")
        if not filename:
            return
        if not filename.lower().endswith(".wrl"):
            filename += ".wrl"
        rw = self.vtk_iren.GetRenderWindow()
        axes_was_on = self._axis_on
        if axes_was_on:
            self._axes_actor.SetVisibility(False)
        rw.Render()
        vrml = vtk.vtkVRMLExporter()
        vrml.SetRenderWindow(rw)
        vrml.SetFileName(filename)
        vrml.SetSpeed(5.5)
        vrml.Write()
        if axes_was_on:
            self._axes_actor.SetVisibility(True)
        rw.Render()

    def on_actionExit(self):
        QApplication.instance().quit()

    def import_shape_from_disc(self, filestr):
        self.shp.read_trunc(filestr)
        self.shp.center_to_zero()
        self.cscale = 1.0 / self._cscale_calc()
        self._axes_actor.SetTotalLength(self._scale_fac / self.cscale,
                                        self._scale_fac / self.cscale,
                                        self._scale_fac / self.cscale)
        self.shp.set_needs_updating(True)
        if self.shp.L_max() < 6:
            self.shp.set_L_max(6)
        self._synchronize_shape2sliders()
        self._update_comboBox_sf()
        self.comboBox_sf.blockSignals(True)
        self.comboBox_sf.setCurrentIndex(0)
        self.comboBox_sf.blockSignals(False)
        self.shp.zero_scalar_field()
        self._update_vtkwindow_full()
        self._ren1.ResetCamera()
        self.vtk_iren.GetRenderWindow().Render()

    # ---- drag & drop ------------------------------------------------------
    def _first_local_file(self, mime):
        if mime.hasUrls():
            for url in mime.urls():
                p = url.toLocalFile()
                if p:
                    return p
        return None

    def dragEnterEvent(self, event):
        if self._first_local_file(event.mimeData()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        p = self._first_local_file(event.mimeData())
        if p:
            print("Opening:", p)
            self.import_shape_from_disc(p)

    def closeEvent(self, event):
        # ensure the VTK interactor/render window shut down cleanly
        try:
            self.vtk_iren.GetRenderWindow().Finalize()
            self.vtk_iren.TerminateApp()
        except Exception:
            pass
        super().closeEvent(event)


def main():
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv)
    win = ShapeWindow()
    win.show()
    # allow loading a .shp3 passed on the command line / file association
    if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
        win.import_shape_from_disc(sys.argv[1])
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
