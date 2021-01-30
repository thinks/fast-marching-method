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

#ifndef FAST_MARCHING_METHOD_TEST_UTIL_HPP_INCLUDED
#define FAST_MARCHING_METHOD_TEST_UTIL_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


namespace util {

template<std::size_t N>
class IndexIterator;

template<typename S, std::size_t N>
struct ScalarDimensionPair
{
  typedef S ScalarType;
  static constexpr std::size_t kDimension = N;
};


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


//! Returns an array with all elements initialized to have value @a v.
template<std::size_t N, typename T> inline
std::array<T, N> FilledArray(T const v)
{
  using namespace std;

  auto r = array<T, N>{};
  fill(begin(r), end(r), v);
  return r;
}


//! Returns the product of the elements in array @a a.
//! Note: Not checking for integer overflow here!
template<std::size_t N> inline
std::size_t LinearSize(std::array<std::size_t, N> const& a)
{
  using namespace std;

  return accumulate(begin(a), end(a), size_t{1}, multiplies<size_t>());
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


//! Returns a pair where the first element is true if the provided function
//! (@a func) threw an exception of the expected type. The second element is
//! set to the exception message.
//!
//! Exceptions that are not of the expected type are re-thrown.
//!
//! E - The expected exception type.
//! F - A callable type.
template<typename E, typename F>
std::pair<bool, std::string> FunctionThrows(F const func)
{
  using namespace std;

  auto thrown_expected = false;
  auto reason = string();
  try {
    func();
  }
  catch (exception& ex) {
    auto const typed_ex = dynamic_cast<E*>(&ex);
    if (typed_ex != nullptr) {
      thrown_expected = true;
      reason = typed_ex->what();
    }
    else {
      throw; // Unexpected exception type, re-throw.
    }
  }

  return {thrown_expected, reason};
}


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


//! Returns true if @a index is inside @a size, otherwise false.
template<std::size_t N>
bool Inside(
  std::array<std::int32_t, N> const& index,
  std::array<std::size_t, N> const& size)
{
  using namespace std;

  static_assert(N > 0, "invalid dimensionality");

  for (auto i = size_t{0}; i < N; ++i) {
    // Cast is safe since we check that index[i] is greater than or
    // equal to zero first.
    if (!(int32_t{0} <= index[i] && static_cast<size_t>(index[i]) < size[i])) {
      return false;
    }
  }
  return true;
}


//! Returns a list of face neighbor offset for an N-dimensional cell.
//! In 2D this is the 4-neighborhood, in 3D the 6-neighborhood, etc.
template<std::size_t N>
std::array<std::array<std::int32_t, N>, 2 * N> FaceNeighborOffsets()
{
  using namespace std;

  static_assert(N > 0, "dimensionality cannot be zero");

  auto offsets = array<array<int32_t, N>, size_t{2} * N>{};
  for (auto i = size_t{0}; i < N; ++i) {
    for (auto j = size_t{0}; j < N; ++j) {
      if (j == i) {
        offsets[2 * i + 0][j] = int32_t{+1};
        offsets[2 * i + 1][j] = int32_t{-1};
      }
      else {
        offsets[2 * i + 0][j] = int32_t{0};
        offsets[2 * i + 1][j] = int32_t{0};
      }
    }
  }
  return offsets;
}


//! DOCS
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
template<typename T, std::size_t N, typename D>
void HyperSphereBoundaryCells(
  std::array<T, N> const& center,
  T const radius,
  std::array<std::size_t, N> const& grid_size,
  std::array<T, N> const& grid_spacing,
  D const distance_modifier,
  std::size_t const dilation_pass_count,
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

  auto foreground_indices = vector<array<int32_t, N>>{};

  // Visit each cell exactly once.
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
      foreground_indices.push_back(index);
    }

    // Update ground truth for all cells.
    if (distance_ground_truth_grid != nullptr) {
      distance_ground_truth_grid->Cell(index) = cell_distance;
    }

    index_iter.Next();
  }

  if (dilation_pass_count > 0 && !foreground_indices.empty()) {
    enum class LabelCell : uint8_t {
      kBackground = uint8_t{0},
      kForeground
    };

    auto label_buffer =
      vector<LabelCell>(LinearSize(grid_size), LabelCell::kBackground);
    auto label_grid = Grid<LabelCell, N>(grid_size, label_buffer.front());
    for (auto const& foreground_index : foreground_indices) {
      label_grid.Cell(foreground_index) = LabelCell::kForeground;
    }

    auto const face_neighbor_offsets = util::FaceNeighborOffsets<N>();
    auto const neighbor_offset_begin = begin(face_neighbor_offsets);
    auto const neighbor_offset_end = end(face_neighbor_offsets);

    auto old_foreground_indices = foreground_indices;
    auto new_foreground_indices = vector<array<int32_t, N>>{};
    for (auto i = size_t{0}; i < dilation_pass_count; ++i) {
      for (auto const foreground_index : old_foreground_indices) {
        for (auto neighbor_offset_iter = neighbor_offset_begin;
             neighbor_offset_iter != neighbor_offset_end;
             ++neighbor_offset_iter) {
          auto neighbor_index = foreground_index;
          for (auto i = size_t{0}; i < N; ++i) {
            neighbor_index[i] += (*neighbor_offset_iter)[i];
          }

          if (Inside(neighbor_index, label_grid.size()) &&
              label_grid.Cell(neighbor_index) == LabelCell::kBackground) {
            auto const cell_center = CellCenter(neighbor_index, grid_spacing);
            auto const cell_distance =
              distance_modifier(Distance(center, cell_center) - radius);
            boundary_indices->push_back(neighbor_index);
            boundary_distances->push_back(cell_distance);

            label_grid.Cell(neighbor_index) = LabelCell::kForeground;
            new_foreground_indices.push_back(neighbor_index);
          }
        }
      }
    }

    old_foreground_indices = new_foreground_indices;
    new_foreground_indices = vector<array<int32_t, N>>{};
  }
}


//! DOCS
template<typename T, std::size_t N>
void BoxBoundaryCells(
  std::array<std::int32_t, N> const& box_corner,
  std::array<std::size_t, N> const& box_size,
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>>* boundary_indices,
  std::vector<T>* boundary_distances)
{
  using namespace std;

  auto box_min = box_corner;
  auto box_max = box_corner;
  for (auto i = size_t{0}; i < N; ++i) {
    box_max[i] += box_size[i];
  }

  auto inside =
    [](auto const& index, auto const& box_min, auto const& box_max) {
      for (auto i = size_t{0}; i < N; ++i) {
        if (!(box_min[i] <= index[i] && index[i] <= box_max[i])) {
          return false;
        }
      }
      return true;
    };

  auto index_iter = util::IndexIterator<N>(grid_size);
  while (index_iter.has_next()) {
    auto const index = index_iter.index();
    if (inside(index, box_min, box_max)) {
      for (auto i = size_t{0}; i < N; ++i) {
        if (index[i] == box_min[i] || index[i] == box_max[i]) {
          boundary_indices->push_back(index);
          break;
        }
      }
    }
    index_iter.Next();
  }

  *boundary_distances = vector<T>(boundary_indices->size(), T{0});
}

} // namespace util

#endif // FAST_MARCHING_METHOD_TEST_UTIL_HPP_INCLUDED
