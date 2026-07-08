"""Non-interactive smoke test: build the window, exercise it, save screenshots."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "app"))
sys.path.insert(0, ROOT)

import shape_app  # imports shp_core first (CRT ordering), then PyQt5 + vtk
import vtk
from PyQt5.QtWidgets import QApplication

OUT = os.path.join(ROOT, "scripts", "out")
os.makedirs(OUT, exist_ok=True)


def snap(win, name):
    rw = win.vtk_iren.GetRenderWindow()
    rw.Render()
    w2i = vtk.vtkWindowToImageFilter()
    w2i.SetInput(rw)
    w2i.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(os.path.join(OUT, name))
    writer.SetInputConnection(w2i.GetOutputPort())
    writer.Write()
    print("saved", name)


def main():
    app = QApplication(sys.argv)
    win = shape_app.ShapeWindow()
    win.resize(1000, 800)
    win.show()
    for _ in range(5):
        app.processEvents()

    snap(win, "01_sphere.png")

    # deform: push a couple of z coefficients (mimics dragging sliders)
    win.verticalSliderx05.setValue(70)
    win.verticalSliderx09xx.setValue(30)
    for _ in range(5):
        app.processEvents()
    snap(win, "02_deformed.png")

    # colour by analytical mean curvature
    idx = win.comboBox_sf.findText("mean curvature (analytical)")
    if idx >= 0:
        win.comboBox_sf.setCurrentIndex(idx)
    for _ in range(5):
        app.processEvents()
    snap(win, "03_curvature.png")

    # load a test shape (discocyte / red blood cell)
    shp = os.path.join(ROOT, "data", "test_shapes", "RBC_and_Vesicles", "discocyte.shp3")
    if os.path.isfile(shp):
        win.import_shape_from_disc(shp)
        for _ in range(5):
            app.processEvents()
        snap(win, "04_discocyte.png")

    # export sanity: write STL + PLY of current shape
    win.polydata  # ensure exists
    stl = vtk.vtkSTLWriter(); stl.SetInputData(win.polydata)
    stl.SetFileName(os.path.join(OUT, "discocyte.stl")); stl.Write()
    print("wrote STL")

    print("SMOKE OK")
    win.close()
    app.quit()


if __name__ == "__main__":
    main()
