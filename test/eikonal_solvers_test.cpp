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

#include "../include/thinks/fast_marching_method/fast_marching_method.hpp"
#include "./util.hpp"

namespace {

// Fixtures.

template <typename T>
class UniformSpeedEikonalSolverTest : public ::testing::Test {
 protected:
  virtual ~UniformSpeedEikonalSolverTest() {}
};

template <typename T>
class HighAccuracyUniformSpeedEikonalSolverTest : public ::testing::Test {
 protected:
  virtual ~HighAccuracyUniformSpeedEikonalSolverTest() {}
};

template <typename T>
class VaryingSpeedEikonalSolverTest : public ::testing::Test {
 protected:
  virtual ~VaryingSpeedEikonalSolverTest() {}
};

template <typename T>
class HighAccuracyVaryingSpeedEikonalSolverTest : public ::testing::Test {
 protected:
  virtual ~HighAccuracyVaryingSpeedEikonalSolverTest() {}
};

template <typename T>
class DistanceSolverTest : public ::testing::Test {
 protected:
  virtual ~DistanceSolverTest() {}
};

// Associate types with fixtures.

typedef ::testing::Types<
    util::ScalarDimensionPair<float, 2>, util::ScalarDimensionPair<float, 3>,
    util::ScalarDimensionPair<float, 4>, util::ScalarDimensionPair<double, 2>,
    util::ScalarDimensionPair<double, 3>, util::ScalarDimensionPair<double, 4>>
    EikonalSolverTypes;

TYPED_TEST_SUITE(UniformSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_SUITE(HighAccuracyUniformSpeedEikonalSolverTest,
                 EikonalSolverTypes);
TYPED_TEST_SUITE(VaryingSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_SUITE(HighAccuracyVaryingSpeedEikonalSolverTest,
                 EikonalSolverTypes);
TYPED_TEST_SUITE(DistanceSolverTest, EikonalSolverTypes);

// UniformSpeedEikonalSolverTest fixture.

TYPED_TEST(UniformSpeedEikonalSolverTest, InvalidGridSpacingThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_grid_spacing_element :
       invalid_grid_spacing_elements) {
    for (auto i = std::size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element;  // Invalid i'th element.
      auto const speed = ScalarType{1};

      auto expected_reason = std::stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
        [[maybe_unused]] auto const eikonal_solver =
            EikonalSolverType(grid_spacing, speed);
        // (void)eikonal_solver;  // pre-C++11
      });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(UniformSpeedEikonalSolverTest, InvalidSpeedThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = invalid_speed;  // Invalid speed!

    auto expected_reason = std::stringstream();
    expected_reason << "invalid speed: " << speed;

    // Act.
    auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
      [[maybe_unused]] auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed);
    });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

// HighAccuracyUniformSpeedEikonalSolverTest fixture.

TYPED_TEST(HighAccuracyUniformSpeedEikonalSolverTest,
           InvalidGridSpacingThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyUniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_grid_spacing_element :
       invalid_grid_spacing_elements) {
    for (auto i = std::size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element;  // Invalid i'th element.
      auto const speed = ScalarType{1};

      auto expected_reason = std::stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
        [[maybe_unused]] auto const eikonal_solver =
            EikonalSolverType(grid_spacing, speed);
      });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(HighAccuracyUniformSpeedEikonalSolverTest, InvalidSpeedThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyUniformSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = invalid_speed;  // Invalid speed!

    auto expected_reason = std::stringstream();
    expected_reason << "invalid speed: " << speed;

    // Act.
    auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
      [[maybe_unused]] auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed);
    });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

// VaryingSpeedEikonalSolverTest fixture.

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidGridSpacingThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_grid_spacing_element :
       invalid_grid_spacing_elements) {
    for (auto i = std::size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element;  // Invalid i'th element.
      auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
      auto const speed_buffer = std::vector<ScalarType>(
          util::LinearSize(speed_grid_size), ScalarType{1});

      auto expected_reason = std::stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<std::invalid_argument>([&]() {
        [[maybe_unused]] auto const eikonal_solver =
            EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

    // Invalid speed in the middle of the buffer.
    auto speed_buffer = std::vector<ScalarType>(
        util::LinearSize(speed_grid_size), ScalarType{1});
    speed_buffer[speed_buffer.size() / 2] = invalid_speed;

    auto expected_reason = std::stringstream();
    expected_reason << "invalid speed: " << invalid_speed;

    // Act.
    auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
      [[maybe_unused]] auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
    });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedBufferThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

  // Buffer size is not linear grid size!
  auto const speed_buffer = std::vector<ScalarType>(
      util::LinearSize(speed_grid_size) - 1, ScalarType{1});

  auto expected_reason = std::stringstream();
  expected_reason << "grid size " << util::ToString(speed_grid_size)
                  << " does not match cell buffer size " << speed_buffer.size();

  // Act.
  auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
    [[maybe_unused]] auto const eikonal_solver =
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedGridSizeThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  for (auto i = std::size_t{0}; i < kDimension; ++i) {
    speed_grid_size[i] = 0;  // Invalid i'th element.
    auto const speed_buffer = std::vector<ScalarType>(
        util::LinearSize(speed_grid_size), ScalarType{1});

    auto expected_reason = std::stringstream();
    expected_reason << "invalid size: " << util::ToString(speed_grid_size);

    // Act.
    auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
      [[maybe_unused]] auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
    });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, IndexOutsideSpeedGridThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  // Speed grid smaller than distance grid!
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{9});
  auto const speed_buffer =
      std::vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

  auto boundary_indices = std::vector<std::array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  auto const boundary_distances = std::vector<ScalarType>(1, ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
    auto const unsigned_distance = fmm::SignedArrivalTime(
        grid_size, boundary_indices, boundary_distances,
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer));
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("speed index outside grid - index:", ft.second.substr(0, 33));
}

// HighAccuracyVaryingSpeedEikonalSolverTest fixture.

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest,
           InvalidGridSpacingThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_grid_spacing_element :
       invalid_grid_spacing_elements) {
    for (auto i = std::size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element;  // Invalid i'th element.
      auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
      auto const speed_buffer = std::vector<ScalarType>(
          util::LinearSize(speed_grid_size), ScalarType{1});

      auto expected_reason = std::stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<std::invalid_argument>([&]() {
        [[maybe_unused]] auto const eikonal_solver =
            EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidSpeedThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

    // Invalid speed in the middle of the buffer.
    auto speed_buffer = std::vector<ScalarType>(
        util::LinearSize(speed_grid_size), ScalarType{1});
    speed_buffer[speed_buffer.size() / 2] = invalid_speed;

    auto expected_reason = std::stringstream();
    expected_reason << "invalid speed: " << invalid_speed;

    // Act.
    auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
      [[maybe_unused]] auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
    });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest,
           InvalidSpeedBufferThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  // Buffer size is not the same as linear grid size!
  auto const speed_buffer = std::vector<ScalarType>(
      util::LinearSize(speed_grid_size) - 1, ScalarType{1});

  auto expected_reason = std::stringstream();
  expected_reason << "grid size " << util::ToString(speed_grid_size)
                  << " does not match cell buffer size " << speed_buffer.size();

  // Act.
  auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
    [[maybe_unused]] auto const eikonal_solver =
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest,
           InvalidSpeedGridSizeThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  for (auto i = std::size_t{0}; i < kDimension; ++i) {
    speed_grid_size[i] = 0;  // Invalid i'th element.
    auto const speed_buffer = std::vector<ScalarType>(
        util::LinearSize(speed_grid_size), ScalarType{1});

    auto expected_reason = std::stringstream();
    expected_reason << "invalid size: " << util::ToString(speed_grid_size);

    // Act.
    auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
      [[maybe_unused]] auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
    });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest,
           IndexOutsideSpeedGridThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
      EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  // Speed grid smaller than distance grid!
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{9});
  auto const speed_buffer =
      std::vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

  auto boundary_indices = std::vector<std::array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  auto const boundary_distances = std::vector<ScalarType>(1, ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
    auto const unsigned_distance = fmm::SignedArrivalTime(
        grid_size, boundary_indices, boundary_distances,
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer));
  });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("speed index outside grid - index:", ft.second.substr(0, 33));
}

// DistanceSolveTest fixture.

TYPED_TEST(DistanceSolverTest, InvalidGridSpacingThrows) {
  typedef typename TypeParam::ScalarType ScalarType;
  static constexpr std::size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::DistanceSolver<ScalarType, kDimension> EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = std::array<ScalarType, 4>{
      {ScalarType{0}, ScalarType{-1},
       std::numeric_limits<ScalarType>::quiet_NaN(), ScalarType(1e-7)}};

  for (auto const invalid_grid_spacing_element :
       invalid_grid_spacing_elements) {
    for (auto i = std::size_t{0}; i < kDimension; ++i) {
      auto const dx = invalid_grid_spacing_element;  // Invalid i'th element.

      auto expected_reason = std::stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(util::FilledArray<kDimension>(dx));

      // Act.
      auto const ft = util::FunctionThrows<std::invalid_argument>([=]() {
        [[maybe_unused]] auto const eikonal_solver = EikonalSolverType(dx);
      });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

}  // namespace
