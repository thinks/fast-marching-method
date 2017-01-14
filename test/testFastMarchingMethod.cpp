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

#include <gtest/gtest.h>

#include <thinks/fastMarchingMethod.hpp>

// TMP!!!
//#include <iostream>


namespace {
namespace util {

//! Returns @a base ^ @a exponent as a compile-time constant.
constexpr std::size_t static_pow(std::size_t base, std::size_t const exponent)
{
  using namespace std;

  // NOTE: Cannot use loops in constexpr functions in C++11, have to use
  // recursion here.
  return exponent == size_t{0} ?
    size_t{1} :
    base * static_pow(base, exponent - 1);
}


//! Returns the product of the elements in array @a a.
//! Note: Not checking for integer overflow here!
template<std::size_t N> inline
std::size_t LinearSize(std::array<std::size_t, N> const& a)
{
  using namespace std;
  return accumulate(begin(a), end(a), size_t{1}, multiplies<size_t>());
}


//! Returns an array with all elements initialized to have value @a v.
template<std::size_t N, typename T> inline
std::array<T, N> FilledArray(T const v)
{
  using namespace std;

  auto r = array<T, N>{};
  fill(begin(r), end(r), v);
  return r;
}


//! Returns array @a as a string.
template<typename T, std::size_t N>
std::string ToString(std::array<T, N> const& a)
{
  using namespace std;

  auto ss = stringstream();
  ss << "[";
  for (auto i = size_t{0}; i < N; ++i) {
    ss << a[i];
    if (i != (N - 1)) {
      ss << ", ";
    }
  }
  ss << "]";

  return ss.str();
}


//! Returns the distance between positions @a u and @a v.
template<typename T, std::size_t N>
T Distance(std::array<T, N> const& u, std::array<T, N> const& v)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "scalar type must be floating point");

  auto distance_squared = T{0};
  for (auto i = size_t{0}; i < N; ++i) {
    auto const delta = u[i] - v[i];
    distance_squared += delta * delta;
  }
  return sqrt(distance_squared);
}


//! Returns the magnitude of the vector @a v.
template<typename T, std::size_t N>
T Magnitude(std::array<T, N> const& v)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "scalar type must be floating point");

  auto mag_squared = T{0};
  for (auto i = size_t{0}; i < N; ++i) {
    mag_squared += v[i] * v[i];
  }
  return sqrt(mag_squared);
}


//! Returns a pair where the first element is true if the provided function
//! (@a func) threw an exception. The second element is set to the exception
//! message.
template<typename E, typename F>
std::pair<bool, std::string> FunctionThrows(F const func)
{
  using namespace std;

  auto thrown = false;
  auto reason = string();
  try {
    func();
  }
  catch (exception& ex)
  {
    thrown = true;
    auto const typed_ex = dynamic_cast<E*>(&ex);
    if (typed_ex == nullptr) {
      auto ss = stringstream();
      ss << "incorrect exception type: '" << typeid(ex).name() << "' "
         << "(expected '" << typeid(E).name() << "')";
      reason = ss.str();
      //cerr << ex.what() << endl;
    }
    else {
      reason = ex.what();
    }
  }

  return {thrown, reason};
}


//! Iterates over a size in N dimensions.
//!
//! Usage:
//!   auto size = std::array<std::size_t, 2>();
//!   size[0] = 2;
//!   size[1] = 2;
//!   auto index_iter = IndexIterator<N>(size);
//!   while (index_iter.has_next()) {
//!     auto index = index_iter.index();
//!     std::cout << "[" << index[0] << ", " << index[1] << "]" << std::endl;
//!     index_iter.Next();
//!   }
//!
//! Output:
//!   [0, 0]
//!   [1, 0]
//!   [0, 1]
//!   [1, 1]
template<std::size_t N>
class IndexIterator
{
public:
  explicit IndexIterator(std::array<std::size_t, N> const& size)
    : size_(size)
    , index_(FilledArray<N>(int32_t{0})) // Start at origin.
    , has_next_(true)
  {
    using namespace std;

    static_assert(N > 0, "must have at least one dimension");

    for (auto const s : size) {
      if (s < 1) {
        throw runtime_error("zero element in size");
      }
    }
  }

  //! Returns the current index of the iterator.
  std::array<std::int32_t, N> index() const
  {
    return index_;
  }

  //! Returns the size over which the iterator is iterating.
  std::array<std::size_t, N> size() const
  {
    return size_;
  }

  //! Returns true if the iterator can be incremented (i.e. Next() can be
  //! called without throwing), otherwise false.
  bool has_next() const
  {
    return has_next_;
  }

  //! Increments the index of the iterator and returns true if it can be
  //! incremented again, otherwise false.
  //!
  //! Throws std::runtime_error if the iterator cannot be incremented.
  bool Next()
  {
    using namespace std;

    if (!has_next_) {
      throw runtime_error("index iterator cannot be incremented");
    }

    auto i = size_t{0};
    while (i < N) {
      if (static_cast<size_t>(index_[i] + 1) < size_[i]) {
        ++index_[i];
        return true;
      }
      else {
        index_[i++] = 0;
      }
    }

    // Nothing was incremented.
    has_next_ = false;
    return false;
  }

private:
  std::array<std::size_t, N> const size_;
  std::array<std::int32_t, N> index_;
  bool has_next_;
};


template<std::size_t N> inline
std::array<std::array<std::int32_t, N>, static_pow(3, N) - 1>
VertexNeighborOffsets()
{
  using namespace std;

  auto neighbor_offsets = array<array<int32_t, N>, static_pow(3, N) - 1>();
  auto offset_index = size_t{0};
  auto index_size = array<size_t, N>{};
  fill(begin(index_size), end(index_size), size_t{3});
  auto index_iter = IndexIterator<N>(index_size);
  while (index_iter.has_next()) {
    auto offset = index_iter.index();
    for_each(begin(offset), end(offset), [](auto& d) { d -= int32_t{1}; });
    if (!all_of(begin(offset), end(offset), [](auto const i){ return i == 0; })) {
      neighbor_offsets[offset_index++] = offset;
    }
    index_iter.Next();
  }
  assert(offset_index == static_pow(3, N) - 1);

  return neighbor_offsets;
}


//! Access a linear array as if it were an N-dimensional grid.
//!
//! Usage:
//!   auto size = std::array<std::size_t, 2>();
//!   size[0] = 2;
//!   size[1] = 2;
//!   auto cells = std::vector<float>(4);
//!   cells[0] = 0.f;
//!   cells[1] = 1.f;
//!   cells[2] = 2.f;
//!   cells[3] = 3.f;
//!   auto grid = Grid(size, cells.front());
//!   std::cout << grid.cell({{0, 0}}) << std::endl;
//!   std::cout << grid.cell({{1, 0}}) << std::endl;
//!   std::cout << grid.cell({{0, 1}}) << std::endl;
//!   std::cout << grid.cell({{1, 1}}) << std::endl;
//!
//! Output:
//!   0
//!   1
//!   2
//!   3
template<typename T, std::size_t N>
class Grid
{
public:
  typedef T CellType;
  typedef std::array<std::size_t, N> SizeType;
  typedef std::array<std::int32_t, N> IndexType;

  Grid(SizeType const& size, T& cells)
    : size_(size)
    , cells_(&cells)
  {
    using namespace std;

    auto stride = size_t{1};
    for (auto i = size_t{1}; i < N; ++i) {
      stride *= size_[i - 1];
      strides_[i - 1] = stride;
    }
  }

  SizeType size() const
  {
    return size_;
  }

  //! Returns a reference to the cell at @a index. No range checking!
  CellType& Cell(IndexType const& index)
  {
    return cells_[LinearIndex_(index)];
  }

  //! Returns a const reference to the cell at @a index. No range checking!
  CellType const& Cell(IndexType const& index) const
  {
    return cells_[LinearIndex_(index)];
  }

private:
  //! Returns a linear (scalar) index into an array representing an
  //! N-dimensional grid for integer coordinate @a index.
  //! Note that this function does not check for integer overflow!
  std::size_t LinearIndex_(IndexType const& index) const
  {
    using namespace std;

    assert(0 <= index[0] && static_cast<size_t>(index[0]) < size_[0]);
    auto k = static_cast<size_t>(index[0]);
    for (auto i = size_t{1}; i < N; ++i) {
      assert(0 <= index[i] && static_cast<size_t>(index[i]) < size_[i]);
      k += index[i] * strides_[i - 1];
    }
    return k;
  }

  std::array<std::size_t, N> const size_;
  std::array<std::size_t, N - 1> strides_;
  CellType* const cells_;
};


//! Returns the center position of the cell at @a index using @a grid_spacing.
template<typename T, std::size_t N> inline
std::array<T, N> CellCenter(
  std::array<std::int32_t, N> const& index,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  auto cell_center = array<T, N>{};
  for (auto i = size_t{0}; i < N; ++i) {
    cell_center[i] = (index[i] + T(0.5)) * grid_spacing[i];
  }
  return cell_center;
}


//! Returns an array of the corner positions of the cell at @a index using
//! @a grid_spacing. The number of corner positions depends on the
//! dimensionality of the cell.
template<typename T, std::size_t N> inline
std::array<std::array<T, N>, static_pow(2, N)> CellCorners(
  std::array<std::int32_t, N> const& index,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  auto cell_corners = array<array<T, N>, static_pow(2, N)>{};
  for (auto i = size_t{0}; i < static_pow(2, N); ++i) {
    auto const bits = bitset<N>(i);
    for (auto k = size_t{0}; k < N; ++k) {
      cell_corners[i][k] =
        (index[k] + static_cast<int32_t>(bits[k])) * grid_spacing[k];
    }
  }
  return cell_corners;
}


//! DOCS
template<typename T, std::size_t N, typename D> inline
void HyperSphereBoundaryCells(
  std::array<T, N> const& center,
  T const radius,
  std::array<std::size_t, N> const& grid_size,
  std::array<T, N> const& grid_spacing,
  D const distance_modifier,
  std::vector<std::array<std::int32_t, N>>* boundary_indices,
  std::vector<T>* boundary_distances,
  std::vector<T>* distance_ground_truth_buffer = nullptr)
{
  using namespace std;

  auto distance_ground_truth_grid = unique_ptr<Grid<T, N>>();
  if (distance_ground_truth_buffer != nullptr) {
    distance_ground_truth_buffer->resize(LinearSize(grid_size));
    distance_ground_truth_grid.reset(
      new Grid<T, N>(grid_size, distance_ground_truth_buffer->front()));
  }

  auto index_iter = IndexIterator<N>(grid_size);
  while (index_iter.has_next()) {
    auto const index = index_iter.index();
    auto const cell_corners = CellCorners(index, grid_spacing);

    auto inside_count = size_t{0};
    auto outside_count = size_t{0};

    for (auto const& cell_corner : cell_corners) {
      auto const d = Distance(center, cell_corner);
      if (d < radius) {
        ++inside_count;
      }
      else {
        ++outside_count;
      }
    }

    auto const cell_center = CellCenter(index, grid_spacing);
    auto const cell_distance =
      distance_modifier(Distance(center, cell_center) - radius);

    if (inside_count > 0 && outside_count > 0) {
      // The inferface passes through this cell so we freeze it.
      boundary_indices->push_back(index);
      boundary_distances->push_back(cell_distance);
    }

    if (distance_ground_truth_grid != nullptr) {
      distance_ground_truth_grid->Cell(index) = cell_distance;
    }

    index_iter.Next();
  }
}


// Numbers below inspired by the paper "ON THE IMPLEMENTATION OF FAST
// MARCHING METHODS FOR 3D LATTICES" by J. Andreas Bærentzen.
template<std::size_t N>
struct PointSourceAccuracyBounds;

template<>
struct PointSourceAccuracyBounds<1>
{
  static constexpr double max_abs_error() { return double{1e-3}; }
  static constexpr double avg_abs_error() { return double{1e-3}; }
  static constexpr double high_accuracy_max_abs_error() { return double{1e-3}; }
  static constexpr double high_accuracy_avg_abs_error() { return double{1e-3}; }
};

template<>
struct PointSourceAccuracyBounds<2>
{
  static constexpr double max_abs_error() { return double{1.48}; }
  static constexpr double avg_abs_error() { return double{0.89}; }
  static constexpr double high_accuracy_max_abs_error() { return double{0.29}; }
  static constexpr double high_accuracy_avg_abs_error() { return double{0.14}; }
};

template<>
struct PointSourceAccuracyBounds<3>
{
  static constexpr double max_abs_error() { return double{1.51}; }
  static constexpr double avg_abs_error() { return double{0.92}; }
  static constexpr double high_accuracy_max_abs_error() { return double{0.28}; }
  static constexpr double high_accuracy_avg_abs_error() { return double{0.07}; }
};

template<>
struct PointSourceAccuracyBounds<4>
{
  static constexpr double max_abs_error() { return double{1.98}; }
  static constexpr double avg_abs_error() { return double{1.27}; }
  static constexpr double high_accuracy_max_abs_error() { return double{0.28}; }
  static constexpr double high_accuracy_avg_abs_error() { return double{0.06}; }
};

} // namespace util


// Fixtures.

template<typename T>
class UniformSpeedEikonalSolverTest : public ::testing::Test {
protected:
  virtual ~UniformSpeedEikonalSolverTest() {}
};


template<typename T>
class HighAccuracyUniformSpeedEikonalSolverTest : public ::testing::Test {
protected:
  virtual ~HighAccuracyUniformSpeedEikonalSolverTest() {}
};


template<typename T>
class VaryingSpeedEikonalSolverTest : public ::testing::Test {
protected:
  virtual ~VaryingSpeedEikonalSolverTest() {}
};


template<typename T>
class HighAccuracyVaryingSpeedEikonalSolverTest : public ::testing::Test {
protected:
  virtual ~HighAccuracyVaryingSpeedEikonalSolverTest() {}
};


template<typename T>
class DistanceSolverTest : public ::testing::Test {
protected:
  virtual ~DistanceSolverTest() {}
};


template<typename T>
class UnsignedDistanceTest : public ::testing::Test {
protected:
  virtual ~UnsignedDistanceTest() {}
};


template<typename T>
class SignedDistanceTest : public ::testing::Test {
protected:
  virtual ~SignedDistanceTest() {}
};


// Associate types with fixtures.

template<typename S, std::size_t N>
struct ScalarDimensionPair
{
  typedef S ScalarType;
  static constexpr std::size_t kDimension = N;
};

typedef ::testing::Types<
  ScalarDimensionPair<float, 1>,
  ScalarDimensionPair<float, 2>,
  ScalarDimensionPair<float, 3>,
  ScalarDimensionPair<float, 4>,
  ScalarDimensionPair<double, 1>,
  ScalarDimensionPair<double, 2>,
  ScalarDimensionPair<double, 3>,
  ScalarDimensionPair<double, 4>> EikonalSolverTypes;

typedef ::testing::Types<
  ScalarDimensionPair<float, 1>,
  ScalarDimensionPair<float, 2>,
  ScalarDimensionPair<float, 3>,
  ScalarDimensionPair<float, 4>,
  ScalarDimensionPair<double, 1>,
  ScalarDimensionPair<double, 2>,
  ScalarDimensionPair<double, 3>,
  ScalarDimensionPair<double, 4>> UnsignedDistanceTypes;

typedef ::testing::Types<
  ScalarDimensionPair<float, 2>,
  ScalarDimensionPair<float, 3>,
  ScalarDimensionPair<float, 4>,
  ScalarDimensionPair<double, 2>,
  ScalarDimensionPair<double, 3>,
  ScalarDimensionPair<double, 4>> SignedDistanceTypes;

TYPED_TEST_CASE(UniformSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(HighAccuracyUniformSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(VaryingSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(HighAccuracyVaryingSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(DistanceSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(UnsignedDistanceTest, UnsignedDistanceTypes);
TYPED_TEST_CASE(SignedDistanceTest, SignedDistanceTypes);


// UniformSpeedEikonalSolverTest fixture.

TYPED_TEST(UniformSpeedEikonalSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element; // Invalid i'th element.
      auto const speed = ScalarType{1};

      auto expected_reason = stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [=]() {
          auto const eikonal_solver = EikonalSolverType(grid_spacing, speed);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(UniformSpeedEikonalSolverTest, InvalidSpeedThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = invalid_speed; // Invalid speed!

    auto expected_reason = stringstream();
    expected_reason << "invalid speed: " << speed;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver = EikonalSolverType(grid_spacing, speed);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}


// HighAccuracyUniformSpeedEikonalSolverTest fixture.

TYPED_TEST(HighAccuracyUniformSpeedEikonalSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyUniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element; // Invalid i'th element.
      auto const speed = ScalarType{1};

      auto expected_reason = stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [=]() {
          auto const eikonal_solver = EikonalSolverType(grid_spacing, speed);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(HighAccuracyUniformSpeedEikonalSolverTest, InvalidSpeedThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyUniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = invalid_speed; // Invalid speed!

    auto expected_reason = stringstream();
    expected_reason << "invalid speed: " << speed;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver = EikonalSolverType(grid_spacing, speed);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}


// VaryingSpeedEikonalSolverTest fixture.

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element; // Invalid i'th element.
      auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
      auto const speed_buffer =
        vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

      auto expected_reason = stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [&]() {
          auto const eikonal_solver =
            EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

    // Invalid speed in the middle of the buffer.
    auto speed_buffer =
      vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});
    speed_buffer[speed_buffer.size() / 2] = invalid_speed;

    auto expected_reason = stringstream();
    expected_reason << "invalid speed: " << invalid_speed;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedBufferThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

  // Buffer size is not linear grid size!
  auto const speed_buffer =
    vector<ScalarType>(util::LinearSize(speed_grid_size) - 1, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "grid size " << util::ToString(speed_grid_size)
     << " does not match cell buffer size " << speed_buffer.size();

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const eikonal_solver =
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedGridSizeThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  for (auto i = size_t{0}; i < kDimension; ++i) {
    speed_grid_size[i] = 0; // Invalid i'th element.
    auto const speed_buffer =
      vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

    auto expected_reason = stringstream();
    expected_reason << "invalid size: "
                    << util::ToString(speed_grid_size);

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, IndexOutsideSpeedGridThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  // Speed grid smaller than distance grid!
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{9});
  auto const speed_buffer =
    vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  auto const boundary_distances = vector<ScalarType>(1, ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("speed index outside grid - index:" , ft.second.substr(0, 33));
}


// HighAccuracyVaryingSpeedEikonalSolverTest fixture.

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element; // Invalid i'th element.
      auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
      auto const speed_buffer =
        vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

      auto expected_reason = stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [&]() {
          auto const eikonal_solver =
            EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidSpeedThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

    // Invalid speed in the middle of the buffer.
    auto speed_buffer =
      vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});
    speed_buffer[speed_buffer.size() / 2] = invalid_speed;

    auto expected_reason = stringstream();
    expected_reason << "invalid speed: " << invalid_speed;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidSpeedBufferThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  // Buffer size is not the same as linear grid size!
  auto const speed_buffer =
    vector<ScalarType>(util::LinearSize(speed_grid_size) - 1, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "grid size " << util::ToString(speed_grid_size)
     << " does not match cell buffer size " << speed_buffer.size();

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const eikonal_solver =
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidSpeedGridSizeThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  for (auto i = size_t{0}; i < kDimension; ++i) {
    speed_grid_size[i] = 0; // Invalid i'th element.
    auto const speed_buffer =
      vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

    auto expected_reason = stringstream();
    expected_reason << "invalid size: " << util::ToString(speed_grid_size);

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, IndexOutsideSpeedGridThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  // Speed grid smaller than distance grid!
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{9});
  auto const speed_buffer =
    vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  auto const boundary_distances = vector<ScalarType>(1, ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("speed index outside grid - index:" , ft.second.substr(0, 33));
}


// DistanceSolveTest fixture.

TYPED_TEST(DistanceSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::DistanceSolver<ScalarType, kDimension> EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 4>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN(),
    ScalarType(1e-7)
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto const dx = invalid_grid_spacing_element; // Invalid i'th element.

      auto expected_reason = stringstream();
      expected_reason
        << "invalid grid spacing: "
        << util::ToString(util::FilledArray<kDimension>(dx));

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [=]() {
          auto const eikonal_solver = EikonalSolverType(dx);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}


// UnsignedDistance fixture.

TYPED_TEST(UnsignedDistanceTest, ZeroElementInGridSizeThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  for (auto i = size_t{0}; i < kDimension; ++i) {
    auto grid_size = util::FilledArray<kDimension>(size_t{10});
    grid_size[i] = 0; // Zero element in i'th position.
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = ScalarType{1};

    auto boundary_indices = vector<array<int32_t, kDimension>>();
    boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{0}));
    auto const boundary_distances = vector<ScalarType>(1, ScalarType{1});

    auto expected_reason = stringstream();
    expected_reason << "invalid size: " << util::ToString(grid_size);

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const unsigned_distance = fmm::UnsignedDistance(
          grid_size,
          boundary_indices,
          boundary_distances,
          EikonalSolverType(grid_spacing, speed));
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(UnsignedDistanceTest, EmptyBoundaryThrows)
{
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
  auto const boundary_indices = vector<array<int32_t, kDimension>>{}; // Empty.
  auto const boundary_distances = vector<ScalarType>{}; // Empty.

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("empty boundary condition", ft.second);
}

TYPED_TEST(UnsignedDistanceTest, FullGridBoundaryIndicesThrows)
{
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
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("empty narrow band", ft.second);
}

TYPED_TEST(UnsignedDistanceTest, DuplicateBoundaryIndicesThrows)
{
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
  boundary_indices.push_back(index_iter.index()); // Same index!

  auto const boundary_distances = vector<ScalarType>(2, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "duplicate boundary index: "
                  << util::ToString(index_iter.index());

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(UnsignedDistanceTest, BoundaryIndexOutsideGridThrows)
{
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
                  << "index: " << util::ToString(boundary_indices.back()) << ", "
                  << "grid size: " << util::ToString(grid_size);

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(UnsignedDistanceTest, BoundaryIndicesAndDistancesSizeMismatchThrows)
{
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
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("boundary indices/distances size mismatch", ft.second);
}

TYPED_TEST(UnsignedDistanceTest, InvalidBoundaryDistanceThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_boundary_distances = array<ScalarType, 4>{{
    ScalarType{-1},
    numeric_limits<ScalarType>::max(),
    -numeric_limits<ScalarType>::max(),
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};
  for (auto const invalid_boundary_distance : invalid_boundary_distances) {
    auto boundary_indices = vector<array<int32_t, kDimension>>();
    auto index_iter = util::IndexIterator<kDimension>(grid_size);
    boundary_indices.push_back(index_iter.index());

    auto const boundary_distances =
      vector<ScalarType>(1, invalid_boundary_distance); // Invalid!

    auto expected_reason = stringstream();
    expected_reason << "invalid boundary distance: "
                    << invalid_boundary_distance;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const unsigned_distance = fmm::UnsignedDistance(
          grid_size,
          boundary_indices,
          boundary_distances,
          EikonalSolverType(grid_spacing, speed));
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(UnsignedDistanceTest, EikonalSolverFailThrows)
{
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
  auto const ft = util::FunctionThrows<runtime_error>(
    [=]() {
      auto const unsigned_distance = fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("invalid arrival time (distance)", ft.second.substr(size_t{0}, 31));
}

TYPED_TEST(UnsignedDistanceTest, DifferentUniformSpeed)
{
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
  auto const unsigned_distance1 = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed1));

  auto const unsigned_distance2 = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed2));

  // Assert.
  // Check that the distance is halved when the speed is halved.
  // Note that the boundary distance is zero, which also passes this check.
  for (auto i = size_t{0}; i < unsigned_distance1.size(); ++i) {
    auto const d1 = unsigned_distance1[i];
    auto const d2 = unsigned_distance2[i];
    ASSERT_LE(fabs(speed1 * d1 - speed2 * d2), ScalarType(1e-3));
  }
}

TYPED_TEST(UnsignedDistanceTest, VaryingSpeed)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{41});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{20}));

  auto boundary_distances = vector<ScalarType>();
  boundary_distances.push_back(ScalarType{0});

  auto const speed_grid_size = grid_size;
  auto speed_buffer = vector<ScalarType>(util::LinearSize(speed_grid_size));
  auto speed_grid =
    util::Grid<ScalarType, kDimension>(speed_grid_size, speed_buffer.front());
  auto const speed = ScalarType(1);
  auto const mirror_speed = ScalarType(2);
  auto speed_index_iter = util::IndexIterator<kDimension>(speed_grid.size());
  while (speed_index_iter.has_next()) {
    auto const index = speed_index_iter.index();
    if (index[0] < boundary_indices[0][0]) {
      speed_grid.Cell(index) = mirror_speed;
    }
    else {
      speed_grid.Cell(index) = speed;
    }
    speed_index_iter.Next();
  }

  // Act.
  auto unsigned_distance = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer));

  // Assert.
  auto distance_grid =
    util::Grid<ScalarType, kDimension>(grid_size, unsigned_distance.front());

  auto distance_index_iter = util::IndexIterator<kDimension>(grid_size);
  while (distance_index_iter.has_next()) {
    auto const index = distance_index_iter.index();
    auto mid = true;
    for (auto i = size_t{1}; i < kDimension; ++i) {
      if (index[i] != grid_size[i] / 2) {
        mid = false;
        break;
      }
    }
    if (index[0] > boundary_indices[0][0] && mid) {
      auto mirror_index = index;
      mirror_index[0] = 2 * boundary_indices[0][0] - index[0];
      auto const distance = distance_grid.Cell(index);
      auto const mirror_distance = distance_grid.Cell(mirror_index);
      ASSERT_NEAR(distance * speed,
                  mirror_distance * mirror_speed,
                  1e-6);
    }
    distance_index_iter.Next();
  }
}

TYPED_TEST(UnsignedDistanceTest, NonUniformGridSpacing)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{31});
  auto grid_spacing = util::FilledArray<kDimension>(ScalarType{0});
  for (auto i = size_t{0}; i < kDimension; ++i) {
    grid_spacing[i] = ScalarType{1} / (i + 1);
  }
  auto const uniform_speed = ScalarType(1);

  auto boundary_indices = vector<array<int32_t, kDimension>>(
    size_t{1}, util::FilledArray<kDimension>(int32_t{15}));
  auto boundary_distances = vector<ScalarType>(size_t{1}, ScalarType{0});

  // Act.
  auto unsigned_distance = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, uniform_speed));

  // Assert.
  auto distance_grid =
    util::Grid<ScalarType, kDimension>(grid_size, unsigned_distance.front());
  auto center_position = util::FilledArray<kDimension>(ScalarType(0));
  for (auto i = size_t{0}; i < kDimension; ++i) {
    center_position[i] =
      (boundary_indices[0][i] + ScalarType(0.5)) * grid_spacing[i];
  }
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
    auto const abs_error = abs(gt_dist - dist);

    typedef util::PointSourceAccuracyBounds<kDimension> Bounds;
    ASSERT_LE(abs_error, ScalarType(Bounds::max_abs_error()));

    index_iter.Next();
  }
}

TYPED_TEST(UnsignedDistanceTest, BoxBoundary)
{
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

  auto boundary_distances = vector<ScalarType>(
    boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto const unsigned_distance = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed));

  // Assert.
  for (auto const d : unsigned_distance) {
    ASSERT_GE(d, ScalarType{0});
  }
}

TYPED_TEST(UnsignedDistanceTest, Checkerboard)
{
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

  auto boundary_distances = vector<ScalarType>(
    boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto const unsigned_distance = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed));

  // Assert.
  // This test exists to ensure that it is possible to run this input.
}

TYPED_TEST(UnsignedDistanceTest, OverlappingCircles)
{
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
    sphere_center1,
    sphere_radius1,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return fabs(d); },
    &sphere_boundary_indices1,
    &sphere_boundary_distances1);
  auto sphere_boundary_indices2 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances2 = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    sphere_center2,
    sphere_radius2,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return fabs(d); },
    &sphere_boundary_indices2,
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
  auto const unsigned_distance = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, uniform_speed));

  // Assert.
  // This test exists to ensure that it is possible to run this input.
  //ASSERT_TRUE(false); // Should this be valid??
}

TYPED_TEST(UnsignedDistanceTest, CircleInsideCircle)
{
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
    sphere_center1,
    sphere_radius1,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return fabs(d); },
    &sphere_boundary_indices1,
    &sphere_boundary_distances1);
  auto sphere_boundary_indices2 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances2 = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    sphere_center2,
    sphere_radius2,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return fabs(d); },
    &sphere_boundary_indices2,
    &sphere_boundary_distances2);

  auto boundary_indices = sphere_boundary_indices1;
  auto boundary_distances = sphere_boundary_distances1;
  for (auto i = size_t{0}; i < sphere_boundary_indices2.size(); ++i) {
    boundary_indices.push_back(sphere_boundary_indices2[i]);
    boundary_distances.push_back(sphere_boundary_distances2[i]);
  }

  // Act.
  auto const unsigned_distance = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, uniform_speed));

  // Assert.
  // This test exists to ensure that it is possible to run this input.
}

TYPED_TEST(UnsignedDistanceTest, PointSourceHighAccuracy)
{
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
  auto boundary_indices = vector<array<int32_t, kDimension>>(
    size_t{1}, util::FilledArray<kDimension>(int32_t{20}));

  auto boundary_distances = vector<ScalarType>(size_t{1}, ScalarType{0});

  // Compute exact distances in vertex neighborhood for high accuracy
  // fast marching.
  auto high_accuracy_boundary_indices = vector<array<int32_t, kDimension>>();
  high_accuracy_boundary_indices.push_back(
    util::FilledArray<kDimension>(int32_t{20})); // Center.
  auto const vtx_neighbor_offsets = util::VertexNeighborOffsets<kDimension>();
  for (auto const& vtx_neighbor_offset : vtx_neighbor_offsets) {
    auto index = boundary_indices[0];
    for (auto i = size_t{0}; i < kDimension; ++i) {
      index[i] += vtx_neighbor_offset[i];
    }
    high_accuracy_boundary_indices.push_back(index);
  }
  auto high_accuracy_boundary_distances = vector<ScalarType>();
  high_accuracy_boundary_distances.push_back(ScalarType{0}); // Center.
  auto center_position = util::FilledArray<kDimension>(ScalarType(0));
  for (auto i = size_t{0}; i < kDimension; ++i) {
    center_position[i] =
      (high_accuracy_boundary_indices[0][i] + ScalarType(0.5)) * grid_spacing[i];
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
  auto unsigned_distance = fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, uniform_speed));
  auto high_accuracy_unsigned_distance = fmm::UnsignedDistance(
    grid_size,
    high_accuracy_boundary_indices,
    high_accuracy_boundary_distances,
    HighAccuracyEikonalSolverType(grid_spacing, uniform_speed));

  // Compute errors.
  auto distance_grid =
    util::Grid<ScalarType, kDimension>(grid_size, unsigned_distance.front());
  auto high_accuracy_distance_grid = util::Grid<ScalarType, kDimension>(
    grid_size, high_accuracy_unsigned_distance.front());

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


// SignedDistance fixture.

TYPED_TEST(SignedDistanceTest, ZeroElementInGridSizeThrows)
{
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
    grid_size[i] = 0; // Zero element in i'th position.

    auto expected_reason = stringstream();
    expected_reason << "invalid size: " << util::ToString(grid_size);

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const signed_distance = fmm::SignedDistance(
          grid_size,
          boundary_indices,
          boundary_distances,
          EikonalSolverType(grid_spacing, speed));
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(SignedDistanceTest, EmptyBoundaryThrows)
{
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
  auto const boundary_indices = vector<array<int32_t, kDimension>>{}; // Empty.
  auto const boundary_distances = vector<ScalarType>{}; // Empty.

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const signed_distance = fmm::SignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("empty boundary condition", ft.second);
}

TYPED_TEST(SignedDistanceTest, FullGridBoundaryIndicesThrows)
{
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
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const signed_distance = fmm::SignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("empty narrow band", ft.second);
}

TYPED_TEST(SignedDistanceTest, DuplicateBoundaryIndicesThrows)
{
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
  boundary_indices.push_back(index_iter.index()); // Same index!

  auto const boundary_distances = vector<ScalarType>(2, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "duplicate boundary index: "
                  << util::ToString(index_iter.index());

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const signed_distance = fmm::SignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(SignedDistanceTest, BoundaryIndexOutsideGridThrows)
{
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
                  << "index: " << util::ToString(boundary_indices.back()) << ", "
                  << "grid size: " << util::ToString(grid_size);

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const signed_distance = fmm::SignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(SignedDistanceTest, BoundaryIndicesAndDistancesSizeMismatchThrows)
{
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
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const signed_distance = fmm::SignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("boundary indices/distances size mismatch", ft.second);
}

TYPED_TEST(SignedDistanceTest, InvalidBoundaryDistanceThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_boundary_distances = array<ScalarType, 3>{{
    numeric_limits<ScalarType>::max(),
    -numeric_limits<ScalarType>::max(),
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};
  for (auto const invalid_boundary_distance : invalid_boundary_distances) {
    auto boundary_indices = vector<array<int32_t, kDimension>>();
    auto index_iter = util::IndexIterator<kDimension>(grid_size);
    boundary_indices.push_back(index_iter.index());

    auto const boundary_distances =
      vector<ScalarType>(1, invalid_boundary_distance); // Invalid!

    auto expected_reason = stringstream();
    expected_reason << "invalid boundary distance: "
                    << invalid_boundary_distance;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const signed_distance = fmm::SignedDistance(
          grid_size,
          boundary_indices,
          boundary_distances,
          EikonalSolverType(grid_spacing, speed));
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(SignedDistanceTest, EikonalSolverFailThrows)
{
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
  auto const ft = util::FunctionThrows<runtime_error>(
    [=]() {
      auto const signed_distance = fmm::SignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("negative discriminant", ft.second);
}

TYPED_TEST(SignedDistanceTest, DifferentUniformSpeed)
{
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
  auto const signed_distance1 = fmm::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed1));

  auto const signed_distance2 = fmm::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
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

TYPED_TEST(SignedDistanceTest, BoxBoundaryHasOnlyInside)
{
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

  auto boundary_distances = vector<ScalarType>(
    boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto const signed_distance = fmm::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed));

  // Assert.
  // Check that we have negative distance inside the box.
  for (auto const d : signed_distance) {
    ASSERT_LE(d, ScalarType{0});
  }
}

TYPED_TEST(SignedDistanceTest, BoxBoundaryWithOutsideCorner)
{
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

  auto boundary_distances = vector<ScalarType>(
    boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto signed_distance = fmm::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed));

  auto const distance_grid = util::Grid<ScalarType, kDimension>(
    grid_size, signed_distance.front());
  auto const mid_cell = util::FilledArray<kDimension>(int32_t{5});
  auto const corner_cell = util::FilledArray<kDimension>(int32_t{9});

  // Assert.
  // Check that we have negative distance inside the box and positive distance
  // in the corner.
  ASSERT_LT(distance_grid.Cell(mid_cell), ScalarType{0});
  ASSERT_GT(distance_grid.Cell(corner_cell), ScalarType{0});
}

TYPED_TEST(SignedDistanceTest, Checkerboard)
{
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

  auto boundary_distances = vector<ScalarType>(
    boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto signed_distance = fmm::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
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
    }
    else {
      ASSERT_LE(distance_grid.Cell(index), ScalarType{0});
    }

    distance_iter.Next();
  }
}

TYPED_TEST(SignedDistanceTest, OverlappingCircles)
{
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
    sphere_center1,
    sphere_radius1,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return d; },
    &sphere_boundary_indices1,
    &sphere_boundary_distances1);
  auto sphere_boundary_indices2 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances2 = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    sphere_center2,
    sphere_radius2,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return d; },
    &sphere_boundary_indices2,
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
  auto signed_distance = fmm::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
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
    auto index = sphere_center_index1; // Get non-x components.
    index[0] = x; // Set x from loop.
    ASSERT_LE(distance_grid.Cell(index), ScalarType{0});
  }
}

TYPED_TEST(SignedDistanceTest, CircleInsideCircleThrows)
{
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
    sphere_center1,
    sphere_radius1,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return d; },
    &sphere_boundary_indices1,
    &sphere_boundary_distances1);
  auto sphere_boundary_indices2 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances2 = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    sphere_center2,
    sphere_radius2,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return d; },
    &sphere_boundary_indices2,
    &sphere_boundary_distances2);

  auto boundary_indices = sphere_boundary_indices1;
  auto boundary_distances = sphere_boundary_distances1;
  for (auto i = size_t{0}; i < sphere_boundary_indices2.size(); ++i) {
    boundary_indices.push_back(sphere_boundary_indices2[i]);
    boundary_distances.push_back(sphere_boundary_distances2[i]);
  }

  // Act.
  auto const signed_distance = fmm::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, uniform_speed));

  // Assert.
  // This test exists to ensure that it is possible to run this input.

}

TYPED_TEST(SignedDistanceTest, PointSourceHighAccuracy)
{
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
    util::FilledArray<kDimension>(int32_t{20})); // Center.
  auto const vtx_neighbor_offsets = util::VertexNeighborOffsets<kDimension>();
  for (auto const& vtx_neighbor_offset : vtx_neighbor_offsets) {
    auto index = boundary_indices[0];
    for (auto i = size_t{0}; i < kDimension; ++i) {
      index[i] += vtx_neighbor_offset[i];
    }
    high_accuracy_boundary_indices.push_back(index);
  }
  auto high_accuracy_boundary_distances = vector<ScalarType>();
  high_accuracy_boundary_distances.push_back(ScalarType{0}); // Center.
  auto center_position = util::FilledArray<kDimension>(ScalarType(0));
  for (auto i = size_t{0}; i < kDimension; ++i) {
    center_position[i] =
      (high_accuracy_boundary_indices[0][i] + ScalarType(0.5)) * grid_spacing[i];
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
  auto signed_distance = fmm::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, uniform_speed));
  auto high_accuracy_signed_distance = fmm::SignedDistance(
    grid_size,
    high_accuracy_boundary_indices,
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

#if 0 // TMP
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

} // namespace
