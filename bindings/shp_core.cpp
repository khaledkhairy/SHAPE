// pybind11 wrapper around the (unmodified) shape_tools core.
//
// This file is *glue only*: every numeric operation is delegated to the
// original shp_surface / spherical_mesh classes in core/shape_tools.h. The goal
// is byte-for-byte identical math to the original SHAPE C++ application while
// exposing a clean, fast (numpy zero-copy where possible) surface to Python.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "shape_tools.h"

namespace py = pybind11;

class Shape {
public:
    shp_surface* s;

    Shape(int L_max = 24, int dim = 30, int tri_n = 4) {
        s = new shp_surface(L_max, dim, tri_n);
    }
    ~Shape() { delete s; }

    // ---- lifecycle / IO ---------------------------------------------------
    void initialize_to_sphere() { s->initializeToSphere(); }
    void read_trunc(const std::string& p) { s->read_trunc(p.c_str()); }
    void write(const std::string& p) { s->write(p.c_str()); }
    void center_to_zero() { s->center_to_zero(); }

    // ---- basis / mesh reconfiguration ------------------------------------
    void set_L_max(int L) { s->set_L_max(L); }
    void set_new_basis(int L, int dim) { s->set_new_basis(L, dim); }
    void set_new_spherical_mesh(int L, int n) { s->set_new_spherical_mesh(L, n); }

    // ---- introspection ----------------------------------------------------
    int L_max() const { return s->b->L_max; }
    int gdim() const { return s->b->dim; }
    int tri_n() const { return s->sm->tri_n; }
    int n_points() const { return s->sm->n_points; }
    int n_faces() const { return s->sm->n_faces; }
    int n_coeffs() const { return static_cast<int>(s->xc.rows()); }
    bool curv_available() const { return s->sm->_curv_calc; }

    bool get_needs_updating() const { return s->needsUpdating; }
    void set_needs_updating(bool v) { s->needsUpdating = v; }

    // ---- coefficient access (xc/yc/zc are public Eigen column vectors) ----
    double get_xc(int i) const { return s->xc(i); }
    double get_yc(int i) const { return s->yc(i); }
    double get_zc(int i) const { return s->zc(i); }
    void set_xc(int i, double v) { s->xc(i) = v; }
    void set_yc(int i, double v) { s->yc(i) = v; }
    void set_zc(int i, double v) { s->zc(i) = v; }

    // ---- surface synthesis (mirrors shape.cpp update_vtkwindow) -----------
    void update_full() {
        s->needsUpdating = true;
        s->update();
        s->surfaceGen();
        s->update_tri();
        s->needsUpdating = false;
    }
    void update_incremental(int coordix, int cix, double old_val, double new_val) {
        s->needsUpdating = true;
        s->update(coordix, cix, old_val, new_val);
        s->surfaceGen(coordix, cix, old_val, new_val);
        s->update_tri();
        s->needsUpdating = false;
    }

    int self_intersect() { return s->self_intersect(); }
    std::vector<std::string> sf_tags() { return s->sf_tags; }

    // ---- analytical surface properties -----------------------------------
    double A() const { return s->A; }
    double Vol() const { return s->V; }
    double reduced_v() const { return s->v; }
    double Eb() const { return s->Eb; }
    double h() const { return s->h; }
    double T() const { return s->T; }
    // ---- mesh-approximation surface properties ---------------------------
    double smA() const { return s->sm->A; }
    double smV() const { return s->sm->V; }
    double smv() const { return s->sm->v; }
    double smEb() const { return s->sm->Eb; }
    double smh() const { return s->sm->h; }

    // ---- geometry for rendering ------------------------------------------
    py::array_t<double> get_points() {
        const int n = s->sm->n_points;
        py::array_t<double> arr({n, 3});
        auto r = arr.mutable_unchecked<2>();
        for (int i = 0; i < n; i++) {
            r(i, 0) = s->sm->X[i][0];
            r(i, 1) = s->sm->X[i][1];
            r(i, 2) = s->sm->X[i][2];
        }
        return arr;
    }
    // Faces already converted to 0-based indexing (source is 1-based / Matlab).
    py::array_t<long long> get_faces() {
        const int m = s->sm->n_faces;
        py::array_t<long long> arr({m, 3});
        auto r = arr.mutable_unchecked<2>();
        for (int i = 0; i < m; i++) {
            r(i, 0) = s->sm->f0[i] - 1;
            r(i, 1) = s->sm->f1[i] - 1;
            r(i, 2) = s->sm->f2[i] - 1;
        }
        return arr;
    }

    // Scalar field selection, mirroring shape.cpp::update_sf exactly.
    py::array_t<double> compute_scalar_field(int val) {
        const int n = s->sm->n_points;
        if (val <= 0) { for (int i = 0; i < n; i++) s->sm->sf[i] = 1.0; }
        else if (val == 1) { for (int i = 0; i < n; i++) s->sm->sf[i] = s->sm->t(i); }
        else if (val == 2) { for (int i = 0; i < n; i++) s->sm->sf[i] = s->sm->p(i); }
        else if (val == 3) { for (int i = 0; i < n; i++) s->sm->sf[i] = s->sm->Han[i]; }
        else if (val == 4) { for (int i = 0; i < n; i++) s->sm->sf[i] = s->sm->H[i]; }
        else if (val == 5) { for (int i = 0; i < n; i++) s->sm->sf[i] = s->sm->KGan[i]; }
        else if (val == 6) { for (int i = 0; i < n; i++) s->sm->sf[i] = s->sm->X[i][0]; }
        else if (val == 7) { for (int i = 0; i < n; i++) s->sm->sf[i] = s->sm->X[i][1]; }
        else if (val == 8) { for (int i = 0; i < n; i++) s->sm->sf[i] = s->sm->X[i][2]; }
        else { s->sfGen(val - 9); }

        py::array_t<double> arr(n);
        auto r = arr.mutable_unchecked<1>();
        for (int i = 0; i < n; i++) r(i) = s->sm->sf[i];
        return arr;
    }

    void zero_scalar_field() {
        for (int i = 0; i < s->sm->n_points; i++) s->sm->sf[i] = 0.0;
    }
};

PYBIND11_MODULE(shp_core, m) {
    m.doc() = "SHAPE core (SPHARM surface) bound from the original C++ shape_tools";

    py::class_<Shape>(m, "Shape")
        .def(py::init<int, int, int>(),
             py::arg("L_max") = 24, py::arg("dim") = 30, py::arg("tri_n") = 4)
        .def("initialize_to_sphere", &Shape::initialize_to_sphere)
        .def("read_trunc", &Shape::read_trunc, py::arg("path"))
        .def("write", &Shape::write, py::arg("path"))
        .def("center_to_zero", &Shape::center_to_zero)
        .def("set_L_max", &Shape::set_L_max, py::arg("L"))
        .def("set_new_basis", &Shape::set_new_basis, py::arg("L"), py::arg("dim"))
        .def("set_new_spherical_mesh", &Shape::set_new_spherical_mesh, py::arg("L"), py::arg("n"))
        .def("L_max", &Shape::L_max)
        .def("gdim", &Shape::gdim)
        .def("tri_n", &Shape::tri_n)
        .def("n_points", &Shape::n_points)
        .def("n_faces", &Shape::n_faces)
        .def("n_coeffs", &Shape::n_coeffs)
        .def("curv_available", &Shape::curv_available)
        .def("get_needs_updating", &Shape::get_needs_updating)
        .def("set_needs_updating", &Shape::set_needs_updating, py::arg("v"))
        .def("get_xc", &Shape::get_xc, py::arg("i"))
        .def("get_yc", &Shape::get_yc, py::arg("i"))
        .def("get_zc", &Shape::get_zc, py::arg("i"))
        .def("set_xc", &Shape::set_xc, py::arg("i"), py::arg("v"))
        .def("set_yc", &Shape::set_yc, py::arg("i"), py::arg("v"))
        .def("set_zc", &Shape::set_zc, py::arg("i"), py::arg("v"))
        .def("update_full", &Shape::update_full)
        .def("update_incremental", &Shape::update_incremental,
             py::arg("coordix"), py::arg("cix"), py::arg("old_val"), py::arg("new_val"))
        .def("self_intersect", &Shape::self_intersect)
        .def("sf_tags", &Shape::sf_tags)
        .def("A", &Shape::A)
        .def("Vol", &Shape::Vol)
        .def("reduced_v", &Shape::reduced_v)
        .def("Eb", &Shape::Eb)
        .def("h", &Shape::h)
        .def("T", &Shape::T)
        .def("smA", &Shape::smA)
        .def("smV", &Shape::smV)
        .def("smv", &Shape::smv)
        .def("smEb", &Shape::smEb)
        .def("smh", &Shape::smh)
        .def("get_points", &Shape::get_points)
        .def("get_faces", &Shape::get_faces)
        .def("compute_scalar_field", &Shape::compute_scalar_field, py::arg("val"))
        .def("zero_scalar_field", &Shape::zero_scalar_field);
}
