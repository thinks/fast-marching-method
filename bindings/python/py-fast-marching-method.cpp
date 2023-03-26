#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "../../include/thinks/fast_marching_method/fast_marching_method.hpp"

namespace py = pybind11;
namespace fmm = thinks::fast_marching_method;

PYBIND11_MODULE(py_fast_marching_method, m) {
    m.doc() = R"pbdoc(
        Python bindings for the fast marching method
        -----------------------
        .. currentmodule:: py_fast_marching_method
        .. autosummary::
           :toctree: _generate
           add
    )pbdoc";

    m.def("add", [](int i, int j) { return i + j; }, R"pbdoc(
        Add two numbers
        Some other explanation about the function.
    )pbdoc");

    // m.def("signed_arrival_time", &fmm::SignedArrivalTime, R"pbdoc(
    //     TODO
    //     Some other explanation about the function.
    // )pbdoc");
}
