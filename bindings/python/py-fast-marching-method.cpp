// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
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

  m.def(
      "add", [](int i, int j) { return i + j; }, R"pbdoc(
        Add two numbers
        Some other explanation about the function.
    )pbdoc");

  // m.def("signed_arrival_time", &fmm::SignedArrivalTime, R"pbdoc(
  //     Signed arrival time
  //     Some other explanation about the function.
  // )pbdoc");
}
