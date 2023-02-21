// Copyright 2017 Tommy Hinks
//
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

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include "../include/thinks/fastMarchingMethod.hpp"
#include "util.hpp"

// TMP!!!
// #include <iostream>

namespace {

// Fixtures.

template <typename T>
class SignedDistanceTest : public ::testing::Test {
 protected:
  virtual ~SignedDistanceTest() {}
};

// Associate types with fixtures.

template <typename S, std::size_t N>
struct ScalarDimensionPair {
  typedef S ScalarType;
  static constexpr std::size_t kDimension = N;
};

typedef ::testing::Types<
    ScalarDimensionPair<float, 2>, ScalarDimensionPair<float, 3>,
    ScalarDimensionPair<float, 4>, ScalarDimensionPair<double, 2>,
    ScalarDimensionPair<double, 3>, ScalarDimensionPair<double, 4>>
    SignedDistanceTypes;

TYPED_TEST_SUITE(SignedDistanceTest, SignedDistanceTypes, );

// SignedDistance fixture.

TYPED_TEST(SignedDistanceTest, ZeroElementInGridSizeThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};
  auto boundary_indices = vector<array<int32_t, kDimension>>();
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{0}));
  auto const boundary_distances = vector<ScalarType>(1, ScalarType{1});
  for (auto i = size_t{0}; i < kDimension; ++i) {
    auto grid_size = util::FilledArray<kDimension>(size_t{10});
    grid_size[i] = 0;  // Zero element in i'th position.

    auto expected_reason = stringstream();
    expected_reason << "invalid size: " << util::ToString(grid_size);

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>([=]() {
      auto const signed_distance =
          fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                              EikonalSolverType(grid_spacing, speed));
    });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(SignedDistanceTest, EmptyBoundaryThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};
  auto const boundary_indices = vector<array<int32_t, kDimension>>{};  // Empty.
  auto const boundary_distances = vector<ScalarType>{};                // Empty.

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>([=]() {
    auto const signed_distance =
        fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                            EikonalSolverType(grid_spacing, speed));
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("empty boundary condition", ft.second);
}

TYPED_TEST(SignedDistanceTest, FullGridBoundaryIndicesThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  // Add every cell in the grid!
  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  while (index_iter.has_next()) {
    boundary_indices.push_back(index_iter.index());
    index_iter.Next();
  }

  auto const boundary_distances =
      vector<ScalarType>(util::LinearSize(grid_size), ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>([=]() {
    auto const signed_distance =
        fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                            EikonalSolverType(grid_spacing, speed));
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("empty narrow band", ft.second);
}

TYPED_TEST(SignedDistanceTest, DuplicateBoundaryIndicesThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  boundary_indices.push_back(index_iter.index());  // Same index!

  auto const boundary_distances = vector<ScalarType>(2, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "duplicate boundary index: "
                  << util::ToString(index_iter.index());

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>([=]() {
    auto const signed_distance =
        fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                            EikonalSolverType(grid_spacing, speed));
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(SignedDistanceTest, BoundaryIndexOutsideGridThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  // Outside!
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{-1}));

  auto const boundary_distances = vector<ScalarType>(2, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "boundary index outside grid - "
                  << "index: " << util::ToString(boundary_indices.back())
                  << ", "
                  << "grid size: " << util::ToString(grid_size);

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>([=]() {
    auto const signed_distance =
        fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                            EikonalSolverType(grid_spacing, speed));
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(SignedDistanceTest, BoundaryIndicesAndDistancesSizeMismatchThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  index_iter.Next();
  boundary_indices.push_back(index_iter.index());

  // Two indices, three distances.
  auto const boundary_distances = vector<ScalarType>(3, ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>([=]() {
    auto const signed_distance =
        fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                            EikonalSolverType(grid_spacing, speed));
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("boundary indices/distances size mismatch", ft.second);
}

TYPED_TEST(SignedDistanceTest, InvalidBoundaryDistanceThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_boundary_distances = array<ScalarType, 3>{
      {numeric_limits<ScalarType>::max(), -numeric_limits<ScalarType>::max(),
       numeric_limits<ScalarType>::quiet_NaN()}};

  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};
  for (auto const invalid_boundary_distance : invalid_boundary_distances) {
    auto boundary_indices = vector<array<int32_t, kDimension>>();
    auto index_iter = util::IndexIterator<kDimension>(grid_size);
    boundary_indices.push_back(index_iter.index());

    auto const boundary_distances =
        vector<ScalarType>(1, invalid_boundary_distance);  // Invalid!

    auto expected_reason = stringstream();
    expected_reason << "invalid boundary distance: "
                    << invalid_boundary_distance;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>([=]() {
      auto const signed_distance =
          fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                              EikonalSolverType(grid_spacing, speed));
    });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(SignedDistanceTest, EikonalSolverFailThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Eikonal solver cannot fail in one dimension?
  if (kDimension == 1) {
    return;
  }

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  // Create a scenario where solver has a very small value in one direction
  // and a very large value in another. Cannot resolve gradient for
  // this scenario.
  auto boundary_indices = vector<array<int32_t, kDimension>>();
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{0}));
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{1}));
  auto boundary_distances = vector<ScalarType>();
  boundary_distances.push_back(ScalarType{1000});
  boundary_distances.push_back(ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<runtime_error>([=]() {
    auto const signed_distance =
        fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                            EikonalSolverType(grid_spacing, speed));
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("invalid arrival time (distance)", ft.second.substr(size_t{0}, 31));
}

TYPED_TEST(SignedDistanceTest, DifferentUniformSpeed) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{5}));

  auto boundary_distances = vector<ScalarType>();
  boundary_distances.push_back(ScalarType{0});

  auto const speed1 = ScalarType{1};
  auto const speed2 = ScalarType{2};

  // Act.
  auto const signed_distance1 =
      fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                          EikonalSolverType(grid_spacing, speed1));

  auto const signed_distance2 =
      fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                          EikonalSolverType(grid_spacing, speed2));

  // Assert.
  // Check that the distance is halved when the speed is halved.
  // Note that the boundary distance is zero, which also passes this check.
  for (auto i = size_t{0}; i < signed_distance1.size(); ++i) {
    auto const d1 = signed_distance1[i];
    auto const d2 = signed_distance2[i];
    ASSERT_LE(fabs(speed1 * d1 - speed2 * d2), ScalarType(1e-3));
  }
}

TYPED_TEST(SignedDistanceTest, BoxBoundaryHasOnlyInside) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  if (kDimension == 1) {
    return;
  }

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  while (index_iter.has_next()) {
    auto const index = index_iter.index();
    for (auto i = size_t{0}; i < kDimension; ++i) {
      if (index[i] == 0 || index[i] == grid_size[i] - 1) {
        boundary_indices.push_back(index);
        break;
      }
    }
    index_iter.Next();
  }

  auto boundary_distances =
      vector<ScalarType>(boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto const signed_distance =
      fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                          EikonalSolverType(grid_spacing, speed));

  // Assert.
  // Check that we have negative distance inside the box.
  for (auto const d : signed_distance) {
    ASSERT_LE(d, ScalarType{0});
  }
}

TYPED_TEST(SignedDistanceTest, BoxBoundaryWithOutsideCorner) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  if (kDimension == 1) {
    return;
  }

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  while (index_iter.has_next()) {
    auto const index = index_iter.index();
    // Don't add the top corner to the boundary condition.
    if (index != util::FilledArray<kDimension>(int32_t{9})) {
      for (auto i = size_t{0}; i < kDimension; ++i) {
        if (index[i] == 0 || index[i] == grid_size[i] - 1) {
          boundary_indices.push_back(index);
          break;
        }
      }
    }
    index_iter.Next();
  }

  auto boundary_distances =
      vector<ScalarType>(boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto signed_distance =
      fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                          EikonalSolverType(grid_spacing, speed));

  auto const distance_grid =
      util::Grid<ScalarType, kDimension>(grid_size, signed_distance.front());
  auto const mid_cell = util::FilledArray<kDimension>(int32_t{5});
  auto const corner_cell = util::FilledArray<kDimension>(int32_t{9});

  // Assert.
  // Check that we have negative distance inside the box and positive distance
  // in the corner.
  ASSERT_LT(distance_grid.Cell(mid_cell), ScalarType{0});
  ASSERT_GT(distance_grid.Cell(corner_cell), ScalarType{0});
}

TYPED_TEST(SignedDistanceTest, Checkerboard) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  if (kDimension == 1) {
    return;
  }

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  auto even = [](auto const i) { return i % 2 == 0; };
  while (index_iter.has_next()) {
    auto const index = index_iter.index();
    if (all_of(begin(index), end(index), even) ||
        none_of(begin(index), end(index), even)) {
      boundary_indices.push_back(index);
    }
    index_iter.Next();
  }

  auto boundary_distances =
      vector<ScalarType>(boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto signed_distance =
      fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                          EikonalSolverType(grid_spacing, speed));

  // Assert.
  auto distance_grid =
      util::Grid<ScalarType, kDimension>(grid_size, signed_distance.front());
  auto distance_iter = util::IndexIterator<kDimension>(grid_size);
  while (distance_iter.has_next()) {
    auto const index = distance_iter.index();
    auto border = false;
    for (auto i = size_t{0}; i < kDimension; ++i) {
      if (index[i] == 0 || index[i] + 1 == grid_size[i]) {
        border = true;
        break;
      }
    }

    if (border) {
      ASSERT_GE(distance_grid.Cell(index), ScalarType{0});
    } else {
      ASSERT_LE(distance_grid.Cell(index), ScalarType{0});
    }

    distance_iter.Next();
  }
}

TYPED_TEST(SignedDistanceTest, OverlappingCircles) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{50});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.02));
  auto const uniform_speed = ScalarType{1};

  auto sphere_center1 = util::FilledArray<kDimension>(ScalarType{0.5});
  sphere_center1[0] = ScalarType(0.4);
  auto const sphere_radius1 = ScalarType{0.25};
  auto sphere_center2 = util::FilledArray<kDimension>(ScalarType{0.5});
  sphere_center2[0] = ScalarType(0.6);
  auto const sphere_radius2 = ScalarType{0.25};

  auto sphere_boundary_indices1 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances1 = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
      sphere_center1, sphere_radius1, grid_size, grid_spacing,
      [](ScalarType const d) { return d; }, &sphere_boundary_indices1,
      &sphere_boundary_distances1);
  auto sphere_boundary_indices2 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances2 = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
      sphere_center2, sphere_radius2, grid_size, grid_spacing,
      [](ScalarType const d) { return d; }, &sphere_boundary_indices2,
      &sphere_boundary_distances2);

  auto boundary_indices = sphere_boundary_indices1;
  auto boundary_distances = sphere_boundary_distances1;
  for (auto i = size_t{0}; i < sphere_boundary_indices2.size(); ++i) {
    auto const index = sphere_boundary_indices2[i];
    if (find(begin(boundary_indices), end(boundary_indices), index) ==
        end(boundary_indices)) {
      boundary_indices.push_back(index);
      boundary_distances.push_back(sphere_boundary_distances2[i]);
    }
  }

  // Act.
  auto signed_distance =
      fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                          EikonalSolverType(grid_spacing, uniform_speed));

  // Assert.
  auto distance_grid =
      util::Grid<ScalarType, kDimension>(grid_size, signed_distance.front());

  auto sphere_center_index1 = util::FilledArray<kDimension>(int32_t{0});
  auto sphere_center_index2 = util::FilledArray<kDimension>(int32_t{0});
  for (auto i = size_t{0}; i < kDimension; ++i) {
    sphere_center_index1[i] =
        static_cast<int32_t>(sphere_center1[i] / grid_spacing[i]);
    sphere_center_index2[i] =
        static_cast<int32_t>(sphere_center2[i] / grid_spacing[i]);
  }

  auto x_begin = min(sphere_center_index1[0], sphere_center_index2[0]);
  auto x_end = max(sphere_center_index1[0], sphere_center_index2[0]);
  for (auto x = x_begin; x <= x_end; ++x) {
    auto index = sphere_center_index1;  // Get non-x components.
    index[0] = x;                       // Set x from loop.
    ASSERT_LE(distance_grid.Cell(index), ScalarType{0});
  }
}

TYPED_TEST(SignedDistanceTest, CircleInsideCircleThrows) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{50});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.02));
  auto const uniform_speed = ScalarType{1};

  auto const sphere_center1 = util::FilledArray<kDimension>(ScalarType{0.5});
  auto const sphere_radius1 = ScalarType(0.25);
  auto const sphere_center2 = util::FilledArray<kDimension>(ScalarType{0.5});
  auto const sphere_radius2 = ScalarType(0.45);

  auto sphere_boundary_indices1 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances1 = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
      sphere_center1, sphere_radius1, grid_size, grid_spacing,
      [](ScalarType const d) { return d; }, &sphere_boundary_indices1,
      &sphere_boundary_distances1);
  auto sphere_boundary_indices2 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances2 = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
      sphere_center2, sphere_radius2, grid_size, grid_spacing,
      [](ScalarType const d) { return d; }, &sphere_boundary_indices2,
      &sphere_boundary_distances2);

  auto boundary_indices = sphere_boundary_indices1;
  auto boundary_distances = sphere_boundary_distances1;
  for (auto i = size_t{0}; i < sphere_boundary_indices2.size(); ++i) {
    boundary_indices.push_back(sphere_boundary_indices2[i]);
    boundary_distances.push_back(sphere_boundary_distances2[i]);
  }

  // Act.
  auto const signed_distance =
      fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                          EikonalSolverType(grid_spacing, uniform_speed));

  // Assert.
  // This test exists to ensure that it is possible to run this input.
}

TYPED_TEST(SignedDistanceTest, PointSourceHighAccuracy) {
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;
  typedef fmm::HighAccuracyUniformSpeedEikonalSolver<ScalarType, kDimension>
      HighAccuracyEikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{41});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const uniform_speed = ScalarType(1);

  // Simple point boundary for regular fast marching.
  auto boundary_indices = vector<array<int32_t, kDimension>>();
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{20}));

  auto boundary_distances = vector<ScalarType>();
  boundary_distances.push_back(ScalarType{0});

  // Compute exact distances in vertex neighborhood for high accuracy
  // fast marching.
  auto high_accuracy_boundary_indices = vector<array<int32_t, kDimension>>();
  high_accuracy_boundary_indices.push_back(
      util::FilledArray<kDimension>(int32_t{20}));  // Center.
  auto const vtx_neighbor_offsets = util::VertexNeighborOffsets<kDimension>();
  for (auto const& vtx_neighbor_offset : vtx_neighbor_offsets) {
    auto index = boundary_indices[0];
    for (auto i = size_t{0}; i < kDimension; ++i) {
      index[i] += vtx_neighbor_offset[i];
    }
    high_accuracy_boundary_indices.push_back(index);
  }
  auto high_accuracy_boundary_distances = vector<ScalarType>();
  high_accuracy_boundary_distances.push_back(ScalarType{0});  // Center.
  auto center_position = util::FilledArray<kDimension>(ScalarType(0));
  for (auto i = size_t{0}; i < kDimension; ++i) {
    center_position[i] =
        (high_accuracy_boundary_indices[0][i] + ScalarType(0.5)) *
        grid_spacing[i];
  }
  for (auto j = size_t{1}; j < high_accuracy_boundary_indices.size(); ++j) {
    auto const& index = high_accuracy_boundary_indices[j];
    auto position = util::FilledArray<kDimension>(ScalarType(0));
    for (auto i = size_t{0}; i < kDimension; ++i) {
      position[i] = (index[i] + ScalarType(0.5)) * grid_spacing[i];
    }
    auto delta = util::FilledArray<kDimension>(ScalarType(0));
    for (auto i = size_t{0}; i < kDimension; ++i) {
      delta[i] = center_position[i] - position[i];
    }
    high_accuracy_boundary_distances.push_back(util::Magnitude(delta));
  }

  // Act.
  auto signed_distance =
      fmm::SignedDistance(grid_size, boundary_indices, boundary_distances,
                          EikonalSolverType(grid_spacing, uniform_speed));
  auto high_accuracy_signed_distance = fmm::SignedDistance(
      grid_size, high_accuracy_boundary_indices,
      high_accuracy_boundary_distances,
      HighAccuracyEikonalSolverType(grid_spacing, uniform_speed));

  // Compute errors.
  auto distance_grid =
      util::Grid<ScalarType, kDimension>(grid_size, signed_distance.front());
  auto high_accuracy_distance_grid = util::Grid<ScalarType, kDimension>(
      grid_size, high_accuracy_signed_distance.front());

  auto dist_abs_errors = vector<ScalarType>();
  auto ha_dist_abs_errors = vector<ScalarType>();

  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  while (index_iter.has_next()) {
    auto const index = index_iter.index();
    auto position = util::FilledArray<kDimension>(ScalarType(0));
    for (auto i = size_t{0}; i < kDimension; ++i) {
      position[i] = (index[i] + ScalarType(0.5)) * grid_spacing[i];
    }
    auto delta = util::FilledArray<kDimension>(ScalarType(0));
    for (auto i = size_t{0}; i < kDimension; ++i) {
      delta[i] = center_position[i] - position[i];
    }
    auto const gt_dist = util::Magnitude(delta);
    auto const dist = distance_grid.Cell(index);
    auto const ha_dist = high_accuracy_distance_grid.Cell(index);
    auto const dist_abs_error = abs(dist - gt_dist);
    auto const ha_dist_abs_error = abs(ha_dist - gt_dist);
    if (gt_dist <= ScalarType{20}) {
      dist_abs_errors.push_back(dist_abs_error);
      ha_dist_abs_errors.push_back(ha_dist_abs_error);
    }

#if 0  // TMP
    dist_abs_error_grid.Cell(index) = dist_abs_error;
    ha_dist_abs_error_grid.Cell(index) = ha_dist_abs_error;
#endif

    index_iter.Next();
  }

  auto max_abs_error = ScalarType{0};
  auto avg_abs_error = ScalarType{0};
  for (auto const& dist_abs_error : dist_abs_errors) {
    max_abs_error = max(max_abs_error, dist_abs_error);
    avg_abs_error += dist_abs_error;
  }
  avg_abs_error /= dist_abs_errors.size();

  auto high_accuracy_max_abs_error = ScalarType{0};
  auto high_accuracy_avg_abs_error = ScalarType{0};
  for (auto const& ha_dist_abs_error : ha_dist_abs_errors) {
    high_accuracy_max_abs_error =
        max(high_accuracy_max_abs_error, ha_dist_abs_error);
    high_accuracy_avg_abs_error += ha_dist_abs_error;
  }
  high_accuracy_avg_abs_error /= ha_dist_abs_errors.size();

  // Assert.
  typedef util::PointSourceAccuracyBounds<kDimension> Bounds;
  ASSERT_LE(max_abs_error, ScalarType(Bounds::max_abs_error()));
  ASSERT_LE(avg_abs_error, ScalarType(Bounds::avg_abs_error()));
  ASSERT_LE(high_accuracy_max_abs_error,
            ScalarType(Bounds::high_accuracy_max_abs_error()));
  ASSERT_LE(high_accuracy_avg_abs_error,
            ScalarType(Bounds::high_accuracy_avg_abs_error()));
}

}  // namespace
