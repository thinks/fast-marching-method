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

#include <array>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../include/thinks/fast_marching_method/fast_marching_method.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  namespace fmm = thinks::fast_marching_method;

  auto grid_size = std::array<std::size_t, 2>{{16, 16}};
  auto grid_spacing = std::array<float, 2>{{1.f / 16, 1.f / 16}};
  auto uniform_speed = 1.f;

  std::cout << "Grid size: " << grid_size[0] << ", " << grid_size[1]
            << std::endl;
  std::cout << "Grid spacing: " << grid_spacing[0] << ", " << grid_spacing[1]
            << std::endl;
  std::cout << "Uniform speed: " << uniform_speed << std::endl;

  // Select points close to a true circle
  auto circle_boundary_indices = std::vector<std::array<int32_t, 2>>{
      {{5, 3}},   {{6, 3}},   {{7, 3}},  {{8, 3}},  {{9, 3}},   {{10, 3}},
      {{4, 4}},   {{5, 4}},   {{10, 4}}, {{11, 4}}, {{3, 5}},   {{4, 5}},
      {{11, 5}},  {{12, 5}},  {{3, 6}},  {{12, 6}}, {{3, 7}},   {{12, 7}},
      {{3, 8}},   {{12, 8}},  {{3, 9}},  {{12, 9}}, {{3, 10}},  {{4, 10}},
      {{11, 10}}, {{12, 10}}, {{4, 11}}, {{5, 11}}, {{10, 11}}, {{11, 11}},
      {{5, 12}},  {{6, 12}},  {{7, 12}}, {{8, 12}}, {{9, 12}},  {{10, 12}},
  };

  auto num_seeds = circle_boundary_indices.size();

  // Specify distances of such points to the true circle
  auto circle_boundary_distances = std::vector<float>{
      0.0417385f, 0.0164635f,  0.0029808f,  0.0029808f,  0.0164635f,
      0.0417385f, 0.0293592f,  -0.0111773f, -0.0111773f, 0.0293592f,
      0.0417385f, -0.0111773f, -0.0111773f, 0.0417385f,  0.0164635f,
      0.0164635f, 0.0029808f,  0.0029808f,  0.0029808f,  0.0029808f,
      0.0164635f, 0.0164635f,  0.0417385f,  -0.0111773f, -0.0111773f,
      0.0417385f, 0.0293592f,  -0.0111773f, -0.0111773f, 0.0293592f,
      0.0417385f, 0.0164635f,  0.0029808f,  0.0029808f,  0.0164635f,
      0.0417385f};

  // auto circle_boundary_distances = std::vector<float>(num_seeds, 0.f);

  auto nan = std::numeric_limits<float>::quiet_NaN();
  auto initial_distances = std::vector<float>(grid_size[0] * grid_size[1], nan);
  for (std::size_t i = 0; i < num_seeds; ++i) {
    auto idx = circle_boundary_indices[i];
    initial_distances[idx[0] + idx[1]* grid_size[0]] =
        circle_boundary_distances[i];
  }

  {
    std::cout << "Initial distance map (seeds):" << std::endl;
    std::size_t idx = 0;
    for (std::size_t j = 0; j < grid_size[1]; ++j) {
      for (std::size_t i = 0; i < grid_size[0]; ++i) {
        std::cout << std::setw(6) << std::fixed << std::setprecision(4)
                  << initial_distances[idx++] << '\t';
      }
      std::cout << std::endl;
    }
  }

  // Compute distance map
  auto arrival_times = fmm::SignedArrivalTime(
      grid_size, circle_boundary_indices, circle_boundary_distances,
      fmm::UniformSpeedEikonalSolver<float, 2>(grid_spacing, uniform_speed));

  {
    std::cout << "Distance map:" << std::endl;
    std::size_t idx = 0;
    for (std::size_t j = 0; j < grid_size[1]; ++j) {
      for (std::size_t i = 0; i < grid_size[0]; ++i) {
        std::cout << std::setw(6) << std::fixed << std::setprecision(4)
                  << arrival_times[idx++] << '\t';
      }
      std::cout << std::endl;
    }
  }
}
