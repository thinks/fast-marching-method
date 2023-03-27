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

template <typename T, std::size_t N>
// std::vector<T>
py::array_t<T> UniformSpeedEikonalSignedArrivalTime(
    std::array<std::size_t, N> const& grid_size,
    std::vector<std::array<std::int32_t, N>> const& boundary_indices,
    std::vector<T> const& boundary_times, std::array<T, N> const& grid_spacing,
    T const uniform_speed) {
  // auto eikonal_solver = fmm::UniformSpeedEikonalSolver<T, N>(grid_spacing,
  // uniform_speed); auto eikonal_solver =
  // fmm::HighAccuracyUniformSpeedEikonalSolver<T, N>(grid_spacing,
  // uniform_speed);
  auto eikonal_solver = fmm::DistanceSolver<T, N>(grid_spacing[0]);
  std::vector<T> arrival_times = fmm::SignedArrivalTime(
      grid_size, boundary_indices, boundary_times, eikonal_solver);
  return py::array_t<T>(grid_size, &arrival_times[0]);
}

PYBIND11_MODULE(py_fast_marching_method, m) {
  m.doc() = R"pbdoc(
        Python bindings for the fast marching method
        -----------------------
        .. currentmodule:: py_fast_marching_method
        .. autosummary::
           :toctree: _generate
           uniform_speed_eikonal_signed_arrival_time
    )pbdoc";

  m.def(
      "add", [](int i, int j) { return i + j; }, R"pbdoc(
        Add two numbers
        Some other explanation about the function.
    )pbdoc");

  m.def("uniform_speed_eikonal_signed_arrival_time",
        &UniformSpeedEikonalSignedArrivalTime<double, 2>, R"pbdoc(
      Signed arrival time
      Some other explanation about the function.
  )pbdoc");
}
