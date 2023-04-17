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
py::array_t<T> UniformSpeedSignedArrivalTime(
    std::array<std::size_t, N> const& py_grid_size,
    std::vector<std::array<std::int32_t, N>> const& py_boundary_indices,
    std::vector<T> const& boundary_times,
    std::array<T, N> const& py_grid_spacing, T const uniform_speed) {
  std::array<std::size_t, N> grid_size;
  std::reverse_copy(py_grid_size.begin(), py_grid_size.end(),
                    grid_size.begin());

  std::array<T, N> grid_spacing;
  std::reverse_copy(py_grid_spacing.begin(), py_grid_spacing.end(),
                    grid_spacing.begin());

  std::vector<std::array<std::int32_t, N>> boundary_indices;
  boundary_indices.reserve(py_boundary_indices.size());
  for (auto& py_boundary_index : py_boundary_indices) {
    std::array<std::int32_t, N> boundary_index;
    std::reverse_copy(py_boundary_index.begin(), py_boundary_index.end(),
                      boundary_index.begin());
    boundary_indices.push_back(boundary_index);
  }

  // auto eikonal_solver = fmm::UniformSpeedEikonalSolver<T, N>
  //   (grid_spacing, uniform_speed);

  auto eikonal_solver = fmm::HighAccuracyUniformSpeedEikonalSolver<T, N>(
      grid_spacing, uniform_speed);

  // auto eikonal_solver = fmm::DistanceSolver<T, N>(grid_spacing[0]);

  std::vector<T> arrival_times = fmm::SignedArrivalTime(
      grid_size, boundary_indices, boundary_times, eikonal_solver);

  return py::array_t<T>(py_grid_size, &arrival_times[0]);
}

template <typename T, std::size_t N>
py::array_t<T> VaryingSpeedSignedArrivalTime(
    std::array<std::size_t, N> const& py_grid_size,
    std::vector<std::array<std::int32_t, N>> const& py_boundary_indices,
    std::vector<T> const& boundary_times,
    std::array<T, N> const& py_grid_spacing,
    py::array_t<T, py::array::c_style | py::array::forcecast> py_speed_buffer) {
  auto py_speed_buffer_flat = py_speed_buffer.reshape({py_speed_buffer.size()});
  auto speed_buffer = py_speed_buffer_flat.template cast<std::vector<T>>();

  std::array<std::size_t, N> grid_size;
  std::reverse_copy(py_grid_size.begin(), py_grid_size.end(),
                    grid_size.begin());

  std::array<T, N> grid_spacing;
  std::reverse_copy(py_grid_spacing.begin(), py_grid_spacing.end(),
                    grid_spacing.begin());

  std::vector<std::array<std::int32_t, N>> boundary_indices;
  boundary_indices.reserve(py_boundary_indices.size());
  for (auto& py_boundary_index : py_boundary_indices) {
    std::array<std::int32_t, N> boundary_index;
    std::reverse_copy(py_boundary_index.begin(), py_boundary_index.end(),
                      boundary_index.begin());
    boundary_indices.push_back(boundary_index);
  }

  // auto eikonal_solver = fmm::VaryingSpeedEikonalSolver<T, N>
  //   (grid_spacing, grid_size, speed_buffer);

  auto eikonal_solver = fmm::HighAccuracyVaryingSpeedEikonalSolver<T, N>(
      grid_spacing, grid_size, speed_buffer);

  std::vector<T> arrival_times = fmm::SignedArrivalTime(
      grid_size, boundary_indices, boundary_times, eikonal_solver);

  return py::array_t<T>(py_grid_size, &arrival_times[0]);
}

PYBIND11_MODULE(py_fast_marching_method, m) {
  m.doc() = R"pbdoc(
        Python bindings for the fast marching method
        -----------------------
        .. currentmodule:: py_fast_marching_method
        .. autosummary::
           :toctree: _generate
           uniform_speed_signed_arrival_time
           varying_speed_signed_arrival_time
    )pbdoc";

  m.def("uniform_speed_signed_arrival_time",
        &UniformSpeedSignedArrivalTime<double, 2>, R"pbdoc(
      Signed arrival time under uniform speed
      https://github.com/thinks/fast-marching-method#high-accuracy-fast-marching-method
  )pbdoc");

  m.def("uniform_speed_signed_arrival_time",
        &UniformSpeedSignedArrivalTime<double, 3>, R"pbdoc(
      Signed arrival time under uniform speed
      https://github.com/thinks/fast-marching-method#high-accuracy-fast-marching-method
  )pbdoc");

  m.def("varying_speed_signed_arrival_time",
        &VaryingSpeedSignedArrivalTime<double, 2>, R"pbdoc(
      Signed arrival time under varying speed
      https://github.com/thinks/fast-marching-method#high-accuracy-fast-marching-method
  )pbdoc");

  m.def("varying_speed_signed_arrival_time",
        &VaryingSpeedSignedArrivalTime<double, 3>, R"pbdoc(
      Signed arrival time under varying speed
      https://github.com/thinks/fast-marching-method#high-accuracy-fast-marching-method
  )pbdoc");
}
