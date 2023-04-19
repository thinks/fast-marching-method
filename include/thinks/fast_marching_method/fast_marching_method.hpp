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

#ifndef INCLUDE_THINKS_FAST_MARCHING_METHOD_FAST_MARCHING_METHOD_HPP_
#define INCLUDE_THINKS_FAST_MARCHING_METHOD_FAST_MARCHING_METHOD_HPP_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <stack>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace thinks {
namespace fast_marching_method {
namespace detail {

// forward decls
template <std::size_t N>
void ThrowIfZeroElementInSize(std::array<std::size_t, N> const& size);

//! Returns the product of the elements in array @a size.
//! Note: Not checking for integer overflow here!
template <std::size_t N>
std::size_t LinearSize(std::array<std::size_t, N> const& size) {
  ThrowIfZeroElementInSize(size);
  return std::accumulate(std::begin(size), std::end(size), size_t{1},
                         std::multiplies<std::size_t>());
}

//! Returns x * x.
//! Note: Not checking for overflow!
template <typename T>
constexpr T Squared(T const x) {
  return x * x;
}

//! Returns 1 / (x * x).
//! Type T must be floating point.
template <typename T>
constexpr T InverseSquared(T const x) {
  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  return T{1} / Squared(x);
}

//! Returns element-wise 1 / (a * a), i.e.
//! { 1 / a[0] * a[0], 1 / a[1] * a[1], ... }
template <typename T, std::size_t N>
std::array<T, N> InverseSquared(std::array<T, N> const& a) {
  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  auto r = std::array<T, N>();
  std::transform(std::begin(a), std::end(a), std::begin(r),
                 [](T const x) { return InverseSquared(x); });
  return r;
}

//! Returns true if @a index is inside @a size, otherwise false.
template <std::size_t N>
bool Inside(std::array<std::int32_t, N> const& index,
            std::array<std::size_t, N> const& size) {
  static_assert(N > 0, "invalid dimensionality");

  for (auto i = std::size_t{0}; i < N; ++i) {
    // Cast is safe since we check that index[i] is greater than or
    // equal to zero first.
    if (!(int32_t{0} <= index[i] &&
          static_cast<std::size_t>(index[i]) < size[i])) {
      return false;
    }
  }
  return true;
}

//! Returns a string representation of the array @a a.
template <typename T, std::size_t N>
std::string ToString(std::array<T, N> const& a) {
  auto ss = std::stringstream();
  ss << "[";
  for (auto i = std::size_t{0}; i < N; ++i) {
    ss << a[i];
    if (i != (N - 1)) {
      ss << ", ";
    }
  }
  ss << "]";
  return ss.str();
}

//! Returns an array where every element is set to @a value.
template <typename T, std::size_t N>
std::array<T, N> FilledArray(T const value) {
  auto a = std::array<T, N>();
  std::fill(std::begin(a), std::end(a), value);
  return a;
}

//! Returns true if @a dx is valid, otherwise false.
template <typename T>
bool ValidGridSpacingElement(T const dx) {
  static_assert(std::is_floating_point<T>::value,
                "grid spacing element type must be floating point");

  // Fail when dx is NaN.
  return dx >= T(1e-6);
}

//! Returns true if @a grid_spacing is valid, otherwise false.
template <typename T, std::size_t N>
bool ValidGridSpacing(std::array<T, N> const& grid_spacing) {
  // All elements must be larger than or equal to a minimum value.
  // Fails if any element is NaN.
  return std::all_of(std::begin(grid_spacing), std::end(grid_spacing),
                     [](auto const dx) { return ValidGridSpacingElement(dx); });
}

//! Returns true if @a speed is valid, otherwise false.
template <typename T>
bool ValidSpeed(T const speed) {
  static_assert(std::is_floating_point<T>::value,
                "speed type must be floating point");

  return speed >= T(1e-6);
}

//! Throws an std::invalid_argument exception if one or more of the
//! elements in @a size is zero.
template <std::size_t N>
void ThrowIfZeroElementInSize(std::array<std::size_t, N> const& size) {
  if (std::find_if(std::begin(size), std::end(size), [](auto const x) {
        return x == std::size_t{0};
      }) != std::end(size)) {
    auto ss = std::stringstream();
    ss << "invalid size: " << ToString(size);
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::invalid_argument exception if the linear size
//! of @a grid_size is not equal to @a cell_buffer_size.
template <std::size_t N>
void ThrowIfInvalidCellBufferSize(std::array<std::size_t, N> const& grid_size,
                                  std::size_t const cell_buffer_size) {
  if (LinearSize(grid_size) != cell_buffer_size) {
    auto ss = std::stringstream();
    ss << "grid size " << ToString(grid_size)
       << " does not match cell buffer size " << cell_buffer_size;
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::invalid_argument exception if @a grid_spacing is invalid.
template <typename T, std::size_t N>
void ThrowIfInvalidGridSpacing(std::array<T, N> const& grid_spacing) {
  if (!ValidGridSpacing(grid_spacing)) {
    auto ss = std::stringstream();
    ss << "invalid grid spacing: " << ToString(grid_spacing);
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::invalid_argument exception if @a speed is invalid.
template <typename T>
void ThrowIfZeroOrNegativeOrNanSpeed(T const speed) {
  if (!ValidSpeed(speed)) {
    auto ss = std::stringstream();
    ss << "invalid speed: " << speed;
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::invalid_argument exception if @a boundary_indices is empty.
template <std::size_t N>
void ThrowIfEmptyBoundaryIndices(
    std::vector<std::array<std::int32_t, N>> const& boundary_indices) {
  if (boundary_indices.empty()) {
    throw std::invalid_argument("empty boundary condition");
  }
}

//! Throws an std::invalid_argument exception if @a boundary_index is not
//! inside @a grid_size.
template <std::size_t N>
void ThrowIfBoundaryIndexOutsideGrid(
    std::array<std::int32_t, N> const& boundary_index,
    std::array<std::size_t, N> const& grid_size) {
  if (!Inside(boundary_index, grid_size)) {
    auto ss = std::stringstream();
    ss << "boundary index outside grid - "
       << "index: " << ToString(boundary_index) << ", "
       << "grid size: " << ToString(grid_size);
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::invalid_argument exception if the flag @a valid is false.
//! @a time is used to construct the exception message.
template <typename T>
void ThrowIfInvalidBoundaryTime(bool const valid, T const time) {
  if (!valid) {
    auto ss = std::stringstream();
    ss << "invalid boundary time: " << time;
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::invalid_argument if the flag @a duplicate is false.
//! @a index is used to construct the exception message.
template <std::size_t N>
void ThrowIfDuplicateBoundaryIndex(bool const duplicate,
                                   std::array<std::int32_t, N> const& index) {
  if (duplicate) {
    auto ss = std::stringstream();
    ss << "duplicate boundary index: " << ToString(index);
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::invalid_argument if the size of @a boundary_indices is
//! not equal to the size of @a boundary_times.
template <typename T, std::size_t N>
void ThrowIfBoundaryIndicesTimesSizeMismatch(
    std::vector<std::array<std::int32_t, N>> const& boundary_indices,
    std::vector<T> const& boundary_times) {
  if (boundary_indices.size() != boundary_times.size()) {
    auto ss = std::stringstream();
    ss << "boundary indices[" << boundary_indices.size() << "] / "
       << "boundary times[" << boundary_times.size() << "] size mismatch";
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::runtime_error exception if @a arrival_time is not valid.
template <typename T, std::size_t N>
void ThrowIfInvalidArrivalTime(T const arrival_time,
                               std::array<int32_t, N> const& index) {
  // Fail when d is NaN.
  if (!(arrival_time >= T{0})) {
    auto ss = std::stringstream();
    ss << "invalid arrival time (distance) " << arrival_time << " at index "
       << ToString(index);
    throw std::runtime_error(ss.str());
  }
}

//! Throws an std::runtime_error exception if @a arrival_time is not valid.
template <typename T, std::size_t N>
void ThrowIfInvalidArrivalTimeWithDetails(
    T const arrival_time, std::array<int32_t, N> const& index, T const speed,
    std::array<T, 3> const& q, T const cell_distance,
    std::array<std::pair<T, std::size_t>, N> const& frozen_neighbor_distances,
    std::size_t const frozen_neighbor_distances_count) {
  // Fail when d is NaN.
  if (!(arrival_time >= T{0})) {
    auto ss = std::stringstream();

    ss << "invalid arrival time (distance) " << arrival_time << " at index "
       << ToString(index) << std::endl;

    // auto const cell_distance = distance_grid.Cell(index);
    for (auto i = std::size_t{0}; i < frozen_neighbor_distances_count; ++i) {
      ss << "[" << i << "]:" << cell_distance << " | "
         << frozen_neighbor_distances[i].first << ", "
         << frozen_neighbor_distances[i].second << std::endl;
    }
    ss << "speed: " << speed << std::endl;
    ss << "q: " << ToString(q) << std::endl;
    ss << "discr: " << q[1] * q[1] - T{4} * q[2] * q[0] << std::endl;
    ss << "int64: "
       << int64_t{-800046} * int64_t{-800046} -
              int64_t{4} * int64_t{1179650} * int64_t{135649}
       << std::endl;

    throw std::runtime_error(ss.str());
  }
}

//! Throws an std::runtime_error exception if @a arrival_time is not valid.
template <typename T, std::size_t N>
void ThrowIfInvalidArrivalTimeWithHighAccuracyDetails(
    T const arrival_time, std::array<int32_t, N> const& index, T const speed,
    std::array<T, 3> const& q, T const cell_distance,
    std::array<std::pair<std::pair<T, T>, std::size_t>, N> const&
        frozen_neighbor_distances,
    std::size_t const frozen_neighbor_distances_count) {
  // Fail when d is NaN.
  if (!(arrival_time >= T{0})) {
    auto ss = std::stringstream();

    ss << "invalid arrival time (distance) " << arrival_time << " at index "
       << ToString(index) << std::endl;

    // auto const cell_distance = distance_grid.Cell(index);
    for (auto i = std::size_t{0}; i < frozen_neighbor_distances_count; ++i) {
      ss << "[" << i << "]:" << cell_distance << " | "
         << frozen_neighbor_distances[i].first.first << ", "
         << frozen_neighbor_distances[i].first.second << std::endl;
    }
    ss << "speed: " << speed << std::endl;
    ss << "q: " << ToString(q) << std::endl;
    ss << "discr: " << q[1] * q[1] - T{4} * q[2] * q[0] << std::endl;
    ss << "int64: "
       << int64_t{-800046} * int64_t{-800046} -
              int64_t{4} * int64_t{1179650} * int64_t{135649}
       << std::endl;

    throw std::runtime_error(ss.str());
  }
}

//! Throws an std::invalid_argument exception if @a speed_index is not
//! inside @a grid_size.
template <std::size_t N>
void ThrowIfSpeedIndexOutsideGrid(
    std::array<std::int32_t, N> const& speed_index,
    std::array<std::size_t, N> const& grid_size) {
  if (!Inside(speed_index, grid_size)) {
    auto ss = std::stringstream();
    ss << "speed index outside grid - "
       << "index: " << ToString(speed_index) << ", "
       << "grid size: " << ToString(grid_size);
    throw std::invalid_argument(ss.str());
  }
}

//! Throws an std::invalid_argument exception if the number of
//! @a boundary_indices is equal to the number of grid cells (given
//! by @a grid_size). Note that we are not checking for duplicates in
//! @a boundary_indices here, we are assuming unique indices. However,
//! duplicate indices are checked elsewhere and will give rise to a separate
//! error.
template <std::size_t N>
void ThrowIfFullGridBoundaryIndices(
    std::vector<std::array<std::int32_t, N>> const& boundary_indices,
    std::array<std::size_t, N> const& grid_size) {
  if (boundary_indices.size() == LinearSize(grid_size)) {
    auto ss = std::stringstream();
    ss << "full grid boundary";
    throw std::invalid_argument(ss.str());
  }
}

//! Returns an std::array that can be used to transform an N-dimensional index
//! into a linear index.
template <std::size_t N>
std::array<std::size_t, N - 1> GridStrides(
    std::array<std::size_t, N> const& grid_size) {
  auto strides = std::array<size_t, N - 1>();
  auto stride = std::size_t{1};
  for (auto i = std::size_t{1}; i < N; ++i) {
    stride *= grid_size[i - 1];
    strides[i - 1] = stride;
  }
  return strides;
}

//! Returns a linear (scalar) index into an std::array representing an
//! N-dimensional grid for integer coordinate @a index.
//! Note: Not checking for integer overflow here!
template <std::size_t N>
std::size_t GridLinearIndex(
    std::array<std::int32_t, N> const& index,
    std::array<std::size_t, N - 1> const& grid_strides) {
  auto k = static_cast<std::size_t>(index[0]);
  for (auto i = std::size_t{1}; i < N; ++i) {
    k += index[i] * grid_strides[i - 1];
  }
  return k;
}

//! Access a linear std::array as if it were an N-dimensional grid.
//! Allows mutating operations on the underlying std::array. The grid does
//! not own the underlying std::array, but is simply an indexing structure.
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
//!   auto grid = Grid<float, 2>(size, cells.front());
//!   std::cout << "cell (0,0): " << grid.cell({{0, 0}}) << std::endl;
//!   std::cout << "cell (1,0): " << grid.cell({{1, 0}}) << std::endl;
//!   std::cout << "cell (0,1): " << grid.cell({{0, 1}}) << std::endl;
//!   std::cout << "cell (1,1): " << grid.cell({{1, 1}}) << std::endl;
//!   std::cout << "-----" << std::endl;
//!   grid.Cell({{0, 1}}) = 5.3;
//!   std::cout << "cell (0,1): " << grid.cell({{0, 1}}) << std::endl;
//!
//! Output:
//!   cell (0,0): 0
//!   cell (1,0): 1
//!   cell (0,1): 2
//!   cell (1,1): 3
//!   -----
//!   cell (0,1): 5.3
template <typename T, std::size_t N>
class Grid {
 public:
  typedef T CellType;
  typedef std::array<std::size_t, N> SizeType;
  typedef std::array<std::int32_t, N> IndexType;

  //! Construct a grid from a given @a size and @a cell_buffer. Does not
  //! take ownership of the cell buffer, it is assumed that this buffer exists
  //! during the life-time of the grid object.
  //!
  //! Preconditions:
  //! - @a cell_buffer is not empty.
  Grid(SizeType const& size, std::vector<T>& cell_buffer)
      : size_(size), strides_(GridStrides(size)), cells_(nullptr) {
    ThrowIfZeroElementInSize(size);
    ThrowIfInvalidCellBufferSize(size, cell_buffer.size());

    assert(!cell_buffer.empty());
    cells_ = &cell_buffer.front();
  }

  //! Returns the size of the grid.
  SizeType size() const { return size_; }

  //! Returns a reference to the cell at @a index. No range checking!
  //!
  //! Preconditions:
  //! - @a index is inside the grid.
  CellType& Cell(IndexType const& index) {
    assert(GridLinearIndex(index, strides_) < LinearSize(size()) &&
           "Precondition");
    return cells_[GridLinearIndex(index, strides_)];
  }

  //! Returns a const reference to the cell at @a index. No range checking!
  //!
  //! Preconditions:
  //! - @a index is inside the grid.
  CellType const& Cell(IndexType const& index) const {
    assert(GridLinearIndex(index, strides_) < LinearSize(size()) &&
           "Precondition");
    return cells_[GridLinearIndex(index, strides_)];
  }

 private:
  std::array<std::size_t, N> const size_;
  std::array<std::size_t, N - 1> const strides_;
  CellType* cells_;
};

//! Access a linear std::array as if it were an N-dimensional grid.
//! Does not allow mutating operations on the underlying std::array. The grid
//! does not own the underlying std::array, but is simply an indexing structure.
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
//!   auto grid = ConstGrid<float, 2>(size, cells.front());
//!   std::cout << "cell (0,0): " << grid.cell({{0, 0}}) << std::endl;
//!   std::cout << "cell (1,0): " << grid.cell({{1, 0}}) << std::endl;
//!   std::cout << "cell (0,1): " << grid.cell({{0, 1}}) << std::endl;
//!   std::cout << "cell (1,1): " << grid.cell({{1, 1}}) << std::endl;
//!
//! Output:
//!   cell (0,0): 0
//!   cell (1,0): 1
//!   cell (0,1): 2
//!   cell (1,1): 3
template <typename T, std::size_t N>
class ConstGrid {
 public:
  typedef T CellType;
  typedef std::array<std::size_t, N> SizeType;
  typedef std::array<std::int32_t, N> IndexType;

  //! Construct a grid from a given @a size and @a cell_buffer. Does not
  //! take ownership of the cell buffer, it is assumed that this buffer exists
  //! during the life-time of the grid object.
  //!
  //! Preconditions:
  //! - @a cell_buffer is not empty.
  ConstGrid(SizeType const& size, std::vector<T> const& cell_buffer)
      : size_(size), strides_(GridStrides(size)), cells_(nullptr) {
    ThrowIfZeroElementInSize(size);
    ThrowIfInvalidCellBufferSize(size, cell_buffer.size());

    assert(!cell_buffer.empty() && "Precondition");
    cells_ = &cell_buffer.front();
  }

  //! Returns the size of the grid.
  SizeType size() const { return size_; }

  //! Returns a const reference to the cell at @a index. No range checking!
  //!
  //! Preconditions:
  //! - @a index is inside the grid.
  CellType const& Cell(IndexType const& index) const {
    assert(GridLinearIndex(index, strides_) < LinearSize(size()) &&
           "Precondition");
    return cells_[GridLinearIndex(index, strides_)];
  }

 private:
  std::array<std::size_t, N> const size_;
  std::array<std::size_t, N - 1> const strides_;
  CellType const* cells_;
};

//! Acceleration structure for keeping track of the smallest distance cell
//! in the narrow band.
template <typename T, std::size_t N>
class NarrowBandStore {
 public:
  typedef T DistanceType;
  typedef std::array<std::int32_t, N> IndexType;
  typedef std::pair<DistanceType, IndexType> ValueType;

  //! Create an empty store.
  NarrowBandStore() {}

  //! Returns true if the store is empty, otherwise false.
  bool empty() const { return min_heap_.empty(); }

  //! Remove the value with the smallest distance from the store and
  //! return it.
  //!
  //! Preconditions:
  //! - The store is not empty (check first with empty()).
  ValueType Pop() {
    assert(!min_heap_.empty() && "Precondition");
    auto const v = min_heap_.top();  // O(1)
    min_heap_.pop();                 // O(log N)
    return v;
  }

  //! Adds @a value to the store.
  void Push(ValueType const& value) {
    min_heap_.push(value);  // O(log N)
  }

 private:
  // Place smaller values at the top of the heap.
  typedef std::priority_queue<ValueType, std::vector<ValueType>,
                              std::greater<ValueType>>
      MinHeap_;

  MinHeap_ min_heap_;
};

//! Returns an std::array of pairs, where each element is the min/max index
//! coordinates in the corresponding dimension.
//!
//! Preconditions:
//! - @a indices is not empty.
template <std::size_t N>
std::array<std::pair<std::int32_t, std::int32_t>, N> BoundingBox(
    std::vector<std::array<std::int32_t, N>> const& indices) {
  static_assert(N > 0, "Dimensionality cannot be zero");

  assert(!indices.empty() && "Precondition");

  // Initialize bounding box in all dimensions.
  auto bbox = std::array<std::pair<int32_t, int32_t>, N>();
  for (auto i = std::size_t{0}; i < N; ++i) {
    bbox[i].first = std::numeric_limits<int32_t>::max();
    bbox[i].second = std::numeric_limits<int32_t>::min();
  }

  // Update with each index in all dimensions.
  for (auto const& index : indices) {
    for (auto i = std::size_t{0}; i < N; ++i) {
      bbox[i].first = std::min(bbox[i].first, index[i]);
      bbox[i].second = std::max(bbox[i].second, index[i]);
    }
  }
  return bbox;
}

//! Returns the hyper volume of the provided N-dimensional
//! bounding box @a bbox. This function does not take grid spacing into
//! account, but rather returns the volume in an index space.
//!
//! Preconditions:
//! - The pairs representing bounds in each dimension store the lower bound
//!   as the first element and the higher bound as the second element.
template <std::size_t N>
std::size_t HyperVolume(
    std::array<std::pair<std::int32_t, std::int32_t>, N> const& bbox) {
  auto hyper_volume = std::size_t{1};
  for (auto i = std::size_t{0}; i < N; ++i) {
    assert(bbox[i].first <= bbox[i].second && "Precondition");
    hyper_volume *= (bbox[i].second - bbox[i].first + 1);
  }
  return hyper_volume;
}

//! Returns true if @a inner_bbox is contained by @a outer_bbox.
//!
//! Preconditions:
//! - The hyper-volume of @a outer_bbox is larger than the hyper-volume of
//!   @a inner_bbox.
//! - The pairs representing bounds in each dimension store the lower bound
//!   as the first element and the higher bound as the second element.
template <std::size_t N>
bool Contains(
    std::array<std::pair<std::int32_t, std::int32_t>, N> const& outer_bbox,
    std::array<std::pair<std::int32_t, std::int32_t>, N> const& inner_bbox) {
  assert(HyperVolume(outer_bbox) >= HyperVolume(inner_bbox) && "Precondition");

  auto contains = true;
  for (auto i = std::size_t{0}; i < N; ++i) {
    assert(outer_bbox[i].first <= outer_bbox[i].second && "Precondition");
    assert(inner_bbox[i].first <= inner_bbox[i].second && "Precondition");
    if (!(outer_bbox[i].first < inner_bbox[i].first &&
          inner_bbox[i].second < outer_bbox[i].second)) {
      contains = false;
      break;
    }
  }
  return contains;
}

//! Returns @a base ^ @a exponent as a compile-time constant.
//! Note: Not checking for integer overflow here!
constexpr std::size_t static_pow(std::size_t const base,
                                 std::size_t const exponent) {
  // Note: Cannot use loops in constexpr functions in C++11, have to use
  // recursion here.
  return exponent == std::size_t{0} ? std::size_t{1}
                                    : base * static_pow(base, exponent - 1);
}

//! Returns a list of face neighbor offset for an N-dimensional cell.
//! In 2D this is the 4-neighborhood, in 3D the 6-neighborhood, etc.
template <std::size_t N>
std::array<std::array<std::int32_t, N>, 2 * N> FaceNeighborOffsets() {
  static_assert(N > 0, "dimensionality cannot be zero");

  auto offsets = std::array<std::array<int32_t, N>, std::size_t{2} * N>{};
  for (auto i = std::size_t{0}; i < N; ++i) {
    for (auto j = std::size_t{0}; j < N; ++j) {
      if (j == i) {
        offsets[2 * i + 0][j] = int32_t{+1};
        offsets[2 * i + 1][j] = int32_t{-1};
      } else {
        offsets[2 * i + 0][j] = int32_t{0};
        offsets[2 * i + 1][j] = int32_t{0};
      }
    }
  }
  return offsets;
}

//! Returns a list of vertex neighbor offset for an N-dimensional cell.
//! In 2D this is the 8-neighborhood, in 3D the 26-neighborhood, etc.
template <std::size_t N>
std::array<std::array<std::int32_t, N>, static_pow(3, N) - 1>
VertexNeighborOffsets() {
  typedef std::array<std::int32_t, N> IndexType;

  static_assert(N > 0, "dimensionality cannot be zero");

  auto offsets = std::array<IndexType, static_pow(3, N) - 1>();
  auto index = IndexType();
  std::fill(std::begin(index), std::end(index), int32_t{0});
  for (auto i = std::size_t{0}; i < offsets.size();) {
    auto offset = index;
    std::for_each(std::begin(offset), std::end(offset),
                  [](auto& d) { d -= int32_t{1}; });
    if (!std::all_of(std::begin(offset), std::end(offset),
                     [](auto const i) { return i == 0; })) {
      offsets[i++] = offset;
    }

    // Next index.
    auto j = std::size_t{0};
    while (j < N) {
      if ((index[j] + 1) < std::int32_t{3}) {
        ++index[j];
        break;
      } else {
        index[j++] = 0;
      }
    }
  }

  return offsets;
}

//! Returns a list of connected components. Each connected component is
//! represented as a non-empty list of indices. The provided
//! @a foreground_indices are used as foreground. The neighborhood used to
//! determine connectivity is given by the two iterators
//! @a neighbor_offset_begin and @a neighbor_offset_end. If @a indices is
//! non-empty there is at least one connected component.
//!
//! Preconditions:
//! - All elements in @a indices are inside @a grid_size.
template <std::size_t N, typename NeighborOffsetIt>
std::vector<std::vector<std::array<std::int32_t, N>>> ConnectedComponents(
    std::vector<std::array<std::int32_t, N>> const& foreground_indices,
    std::array<std::size_t, N> const& grid_size,
    NeighborOffsetIt const neighbor_offset_begin,
    NeighborOffsetIt const neighbor_offset_end) {
  enum class LabelCell : uint8_t {
    kBackground = uint8_t{0},
    kForeground,
    kLabelled
  };

  if (foreground_indices.empty()) {
    return std::vector<std::vector<std::array<int32_t, N>>>();
  }

  auto label_buffer =
      std::vector<LabelCell>(LinearSize(grid_size), LabelCell::kBackground);
  auto label_grid = Grid<LabelCell, N>(grid_size, label_buffer);

  for (auto const& foreground_index : foreground_indices) {
    // Note: We don't check for duplicate indices here since it doesn't
    //       affect the algorithm.
    assert(Inside(foreground_index, label_grid.size()) && "Precondition");
    label_grid.Cell(foreground_index) = LabelCell::kForeground;
  }

  auto connected_components =
      std::vector<std::vector<std::array<std::int32_t, N>>>();
  for (auto const& foreground_index : foreground_indices) {
    assert(Inside(foreground_index, label_grid.size()));
    assert(label_grid.Cell(foreground_index) == LabelCell::kForeground ||
           label_grid.Cell(foreground_index) == LabelCell::kLabelled);

    if (label_grid.Cell(foreground_index) == LabelCell::kForeground) {
      // This index has not already been labelled.
      // Start a new component.
      auto component = std::vector<std::array<int32_t, N>>();
      component.reserve(foreground_indices.size());
      auto neighbor_indices = std::stack<std::array<int32_t, N>>();
      label_grid.Cell(foreground_index) = LabelCell::kLabelled;
      component.push_back(foreground_index);
      neighbor_indices.push(foreground_index);

      // Flood-fill current label.
      while (!neighbor_indices.empty()) {
        auto const top_neighbor_index = neighbor_indices.top();
        neighbor_indices.pop();
        for (auto neighbor_offset_iter = neighbor_offset_begin;
             neighbor_offset_iter != neighbor_offset_end;
             ++neighbor_offset_iter) {
          // Offset neighbor index.
          auto neighbor_index = top_neighbor_index;
          for (auto i = std::size_t{0}; i < N; ++i) {
            neighbor_index[i] += (*neighbor_offset_iter)[i];
          }

          if (Inside(neighbor_index, label_grid.size()) &&
              label_grid.Cell(neighbor_index) == LabelCell::kForeground) {
            // Tag this neighbor as labelled, store in component and add
            // to list of indices whose neighbors we should check.
            label_grid.Cell(neighbor_index) = LabelCell::kLabelled;
            component.push_back(neighbor_index);
            neighbor_indices.push(neighbor_index);
          }
        }
      }

      assert(!component.empty());
      connected_components.push_back(component);
    }
  }

  assert(!connected_components.empty());
  return connected_components;
}

//! Returns @a dilation_grid_index transformed to the distance grid.
template <std::size_t N>
std::array<std::int32_t, N> DistanceGridIndexFromDilationGridIndex(
    std::array<std::int32_t, N> const& dilation_grid_index) {
  auto distance_grid_index = dilation_grid_index;
  std::for_each(std::begin(distance_grid_index), std::end(distance_grid_index),
                [](auto& d) { d -= int32_t{1}; });
  return distance_grid_index;
}

//! Returns @a distance_grid_index transformed to the dilation grid.
template <std::size_t N>
std::array<std::int32_t, N> DilationGridIndexFromDistanceGridIndex(
    std::array<std::int32_t, N> const& distance_grid_index) {
  auto dilation_grid_index = distance_grid_index;
  std::for_each(std::begin(dilation_grid_index), std::end(dilation_grid_index),
                [](auto& d) {
                  d += int32_t{1};
                  assert(d >= int32_t{0});
                });
  return dilation_grid_index;
}

//! Returns a list of dilation bands. A dilation band is defined as a set of
//! cells where each cell has at least one neighbor in @a grid_indices. The
//! neighbor definition is computed using the offsets provided by
//! @a dilation_neighbor_offset_begin and @a dilation_neighbor_offset_end.
//! Furthermore, the cells in a dilation band are connected to each other. In
//! this case the neighborhood is defined by the offsets provided by
//! @a band_neighbor_offset_begin and @a band_neighbor_offset_end.
//!
//! Note that the cells in the dilation bands are defined on a dilation grid
//! that is padded by one in each dimension relative to @a grid_size.
//!
//! Preconditions:
//! - All elements in @a grid_indices are inside @a grid_size.
//! - Neighbor offsets are not larger than one cell in any dimension.
template <std::size_t N, typename DilationNeighborOffsetIt,
          typename BandNeighborOffsetIt>
std::vector<std::vector<std::array<std::int32_t, N>>> DilationBands(
    std::vector<std::array<std::int32_t, N>> const& grid_indices,
    std::array<std::size_t, N> const& grid_size,
    DilationNeighborOffsetIt const dilation_neighbor_offset_begin,
    DilationNeighborOffsetIt const dilation_neighbor_offset_end,
    BandNeighborOffsetIt const band_neighbor_offset_begin,
    BandNeighborOffsetIt const band_neighbor_offset_end) {
  enum class DilationCell : uint8_t {
    kBackground = uint8_t{0},
    kForeground,
    kDilated
  };

  assert(LinearSize(grid_size) > std::size_t{0});

  if (grid_indices.empty()) {
    return std::vector<std::vector<std::array<int32_t, N>>>();
  }

  // Dilation grid is padded one cell in each dimension.
  // Then transform the provided indices to the dilation grid.
  auto dilation_grid_size = grid_size;
  std::for_each(std::begin(dilation_grid_size), std::end(dilation_grid_size),
                [](auto& d) { d += 2; });
  auto dilation_buffer = std::vector<DilationCell>(
      LinearSize(dilation_grid_size), DilationCell::kBackground);
  auto dilation_grid =
      Grid<DilationCell, N>(dilation_grid_size, dilation_buffer);
  auto dilation_grid_indices = std::vector<std::array<int32_t, N>>();
  dilation_grid_indices.reserve(grid_indices.size());
  transform(std::begin(grid_indices), std::end(grid_indices),
            back_inserter(dilation_grid_indices), [=](auto const& grid_index) {
              assert(Inside(grid_index, grid_size));
              auto const dilation_grid_index =
                  DilationGridIndexFromDistanceGridIndex(grid_index);
              assert(Inside(dilation_grid_index, dilation_grid_size));
              return dilation_grid_index;
            });

  // Set foreground from the (transformed) provided indices.
  for (auto const& dilation_grid_index : dilation_grid_indices) {
    // Note: We don't check for duplicate indices here since it doesn't
    //       affect the algorithm.
    assert(Inside(dilation_grid_index, dilation_grid.size()) && "Precondition");
    dilation_grid.Cell(dilation_grid_index) = DilationCell::kForeground;
  }

  // Tag background cells connected to foreground as dilated.
  // We only overwrite background cells here.
  auto dilation_indices = std::vector<std::array<int32_t, N>>();
  dilation_indices.reserve(size_t{2} * dilation_grid_indices.size());
  for (auto const& dilation_grid_index : dilation_grid_indices) {
    assert(dilation_grid.Cell(dilation_grid_index) ==
           DilationCell::kForeground);
    for (auto dilation_neighbor_offset_iter = dilation_neighbor_offset_begin;
         dilation_neighbor_offset_iter != dilation_neighbor_offset_end;
         ++dilation_neighbor_offset_iter) {
      auto neighbor_index = dilation_grid_index;
      for (auto i = std::size_t{0}; i < N; ++i) {
        neighbor_index[i] += (*dilation_neighbor_offset_iter)[i];
      }

      assert(Inside(neighbor_index, dilation_grid.size()) && "Precondition");
      auto& dilation_cell = dilation_grid.Cell(neighbor_index);
      if (dilation_cell == DilationCell::kBackground) {
        dilation_cell = DilationCell::kDilated;
        dilation_indices.push_back(neighbor_index);
      }
    }
  }
  assert(!dilation_indices.empty());

  // Get connected components of dilated cells.
  auto const dilation_bands =
      ConnectedComponents(dilation_indices, dilation_grid_size,
                          band_neighbor_offset_begin, band_neighbor_offset_end);
  assert(!dilation_bands.empty());
  return dilation_bands;
}

//! Since dilation bands are constructed using a vertex neighborhood,
//! not all dilation cells are face-connected to a boundary cell.
//! Given a list of @a dilation_band_indices in dilation grid coordinates
//! (padded by one in each direction), returns a list of narrow band
//! indices in distance grid coordinates. The returned indices are
//! guaranteed to be face-connected to at least one boundary cell
//! in @a boundary_mask_grid.
//!
//! Note that the returned list may be empty. This can happen if all
//! @a dilation_band_indices are on the border of the dilation grid, i.e.
//! outside the distance grid. It also happens if the @a dilation_band_indices
//! list is empty, or if @a boundary_mask_grid has values such that the
//! boundary is not face-connected to any of the dilation indices.
//!
//! It is assumed that the value int8_t{1} is used to tag boundary cells in
//! @a boundary_mask_grid. Also, (transformed) dilation band indices are
//! assumed not to be on a boundary.
template <std::size_t N>
std::vector<std::array<std::int32_t, N>> NarrowBandDilationBandCells(
    std::vector<std::array<std::int32_t, N>> const& dilation_band_indices,
    Grid<uint8_t, N> const& boundary_mask_grid) {
  if (dilation_band_indices.empty()) {
    return std::vector<std::array<int32_t, N>>();
  }

  auto narrow_band_indices = std::vector<std::array<int32_t, N>>();
  narrow_band_indices.reserve(dilation_band_indices.size());
  for (auto const& dilation_grid_index : dilation_band_indices) {
    // Since dilation bands are constructed using a vertex neighborhood,
    // not all dilation cells are face-connected to a boundary cell.
    // We add only those dilation cells that are face-connected to a
    // boundary cell, since this will be required when estimating distance
    // (i.e. solving the eikonal equation).
    auto const distance_grid_index =
        DistanceGridIndexFromDilationGridIndex(dilation_grid_index);

    // If the distance grid index is not inside the boundary mask
    // (i.e. distance) grid it cannot belong to a narrow band.
    if (Inside(distance_grid_index, boundary_mask_grid.size())) {
      assert(boundary_mask_grid.Cell(distance_grid_index) != uint8_t{1});

      // Check for boundary face-neighbors in each dimension.
      // If we find one boundary face-neighbor we are done.
      for (auto i = std::size_t{0}; i < N; ++i) {
        // +1
        auto neighbor_index = distance_grid_index;
        neighbor_index[i] += int32_t{1};
        if (Inside(neighbor_index, boundary_mask_grid.size()) &&
            boundary_mask_grid.Cell(neighbor_index) == uint8_t{1}) {
          narrow_band_indices.push_back(distance_grid_index);
          break;
        }
        // +1 - 2 = -1
        neighbor_index[i] -= int32_t{2};
        if (Inside(neighbor_index, boundary_mask_grid.size()) &&
            boundary_mask_grid.Cell(neighbor_index) == uint8_t{1}) {
          narrow_band_indices.push_back(distance_grid_index);
          break;
        }
      }
    }
  }
  assert(narrow_band_indices.size() <= dilation_band_indices.size());

  return narrow_band_indices;
}

//! Returns a pair of lists:
//! - The first element is the set of cells closest to the boundary that
//!   are on the outside.
//! - The second element is the set of cells closest to the boundary that
//!   are on the inside.
//!
//! One or both lists may be empty. In the case of the outside indices, the
//! returned list may contain duplicates (this is not the case for the inside
//! indices).
//!
//! All returned indices are guaranteed to be inside @a grid_size.
//!
//! Preconditions:
//! - Every element in @a boundary_indices is inside @a grid_size.
template <std::size_t N>
std::pair<std::vector<std::array<std::int32_t, N>>,
          std::vector<std::array<std::int32_t, N>>>
OutsideInsideNarrowBandIndices(
    std::vector<std::array<std::int32_t, N>> const& boundary_indices,
    std::array<std::size_t, N> const& grid_size) {
  auto inside_narrow_band_indices = std::vector<std::array<int32_t, N>>();
  auto outside_narrow_band_indices = std::vector<std::array<int32_t, N>>();
  if (boundary_indices.empty()) {
    return {outside_narrow_band_indices, inside_narrow_band_indices};
  }

  // Compute connected components of boundary cells.
  auto const vtx_neighbor_offsets = VertexNeighborOffsets<N>();
  auto const connected_components = ConnectedComponents(
      boundary_indices, grid_size, begin(vtx_neighbor_offsets),
      std::end(vtx_neighbor_offsets));
  assert(!connected_components.empty());

  // Check if any connected component is contained by another.
  auto const connected_components_size = connected_components.size();
  if (connected_components_size > 1) {
    auto cc_bbox = std::vector<
        std::pair<std::array<std::pair<int32_t, int32_t>, N>, std::size_t>>();
    for (auto const& connected_component : connected_components) {
      auto const bbox = BoundingBox(connected_component);
      cc_bbox.push_back({bbox, HyperVolume(bbox)});
    }
    // Sort by descending area (hyper volume).
    sort(std::begin(cc_bbox), std::end(cc_bbox),
         [](auto const& lhs, auto const& rhs) {
           return lhs.second > rhs.second;
         });
    // A smaller bounding box cannot contain a larger one.
    for (auto i = std::size_t{0}; i < connected_components_size; ++i) {
      auto const& outer_bbox = cc_bbox[i].first;

      // A bounding box that is "flat" in one or more dimensions cannot
      // contain another bounding box.
      auto has_inside = true;
      for (auto k = std::size_t{0}; k < N; ++k) {
        if (outer_bbox[k].first == outer_bbox[k].second) {
          has_inside = false;
          break;
        }
      }

      if (has_inside) {
        for (auto j = i + 1; j < connected_components_size; ++j) {
          auto const& inner_bbox = cc_bbox[j].first;
          if (Contains(outer_bbox, inner_bbox)) {
            throw std::invalid_argument("contained component");
          }
        }
      }
    }
  }

  // Create a mask where:
  // - boundary cells = 1
  // - non-boundary cells = 0
  auto boundary_mask_buffer =
      std::vector<uint8_t>(LinearSize(grid_size), uint8_t{0});
  auto boundary_mask_grid = Grid<uint8_t, N>(grid_size, boundary_mask_buffer);
  for (auto const& boundary_index : boundary_indices) {
    assert(Inside(boundary_index, boundary_mask_grid.size()) && "Precondition");
    boundary_mask_grid.Cell(boundary_index) = uint8_t{1};
  }

  // Check dilation bands of connected boundary components.
  // Dilation bands must be computed per connected component since each
  // component has a separate outer dilation band. If we were to compute
  // dilation bands for all boundary indices at once we would then need to
  // do extra work to figure out if these were outer or inner dilation bands.
  auto const face_neighbor_offsets = FaceNeighborOffsets<N>();
  for (auto const& connected_component : connected_components) {
    auto const dilation_bands = DilationBands(
        connected_component, grid_size, begin(vtx_neighbor_offsets),
        std::end(vtx_neighbor_offsets), begin(face_neighbor_offsets),
        std::end(face_neighbor_offsets));
    assert(!dilation_bands.empty());

    if (dilation_bands.size() == 1) {
      // Only one dilation band means that the connected component has genus
      // zero, i.e. no holes. Thus, the dilation band must define
      // the outside.
      //
      // Note that the outer *dilation band* can never be empty, but the
      // *outer narrow band* can be! The outer narrow band is empty when the
      // whole border of the distance grid is boundary.
      auto const& outer_dilation_band = dilation_bands.front();
      assert(!outer_dilation_band.empty());
      auto const outer_narrow_band_indices =
          NarrowBandDilationBandCells(outer_dilation_band, boundary_mask_grid);
      outside_narrow_band_indices.insert(
          std::end(outside_narrow_band_indices),  // Position.
          begin(outer_narrow_band_indices),
          std::end(outer_narrow_band_indices));
    } else {
      // We have more than one dilation band: one outer and one or more
      // inner. The outer dilation band has the largest bounding box.
      // Note that when we have several dilation bands none of them can be
      // empty. The reasoning is that an empty dilation band requires the
      // whole distance grid to be frozen, in which case there cannot exist
      // an inner area.
      //
      // We compute the bounding boxes in dilation grid coordinates. This is
      // necessary since the entire outer dilation band may not be inside the
      // distance grid.
      auto dilation_band_areas = std::vector<std::pair<size_t, std::size_t>>();
      dilation_band_areas.reserve(dilation_bands.size());
      for (auto i = std::size_t{0}; i < dilation_bands.size(); ++i) {
        [[maybe_unused]] auto const& dilation_band = dilation_bands[i];
        assert(!dilation_band.empty());
        dilation_band_areas.push_back(
            {i, HyperVolume(BoundingBox(dilation_bands[i]))});
      }

      // Sort dilation bands by descending volume. The outer dilation band
      // is then the first element. Note that the outer dilation band area
      // should be strictly larger than the largest inner dilation band area
      // (except in 1D).
      sort(std::begin(dilation_band_areas), std::end(dilation_band_areas),
           [](auto const& lhs, auto const& rhs) {
             return lhs.second > rhs.second;
           });
      assert(N == 1 ||
             dilation_band_areas[0].second > dilation_band_areas[1].second);

      // Outer dilation bands of several connected components may overlap.
      // We are fine with adding an index multiple times to the outside
      // narrow band. The smallest distance will be used first and the rest
      // will be ignored. Worst-case we estimate distances for cells that
      // are not impactful.
      auto const& outer_dilation_band_indices =
          dilation_bands[dilation_band_areas[0].first];
      assert(!outer_dilation_band_indices.empty());
      auto const outer_narrow_band_indices = NarrowBandDilationBandCells(
          outer_dilation_band_indices, boundary_mask_grid);
      assert(none_of(std::begin(outer_narrow_band_indices),
                     std::end(outer_narrow_band_indices),
                     [=](auto const& distance_grid_index) {
                       return !Inside(distance_grid_index, grid_size);
                     }));
      // Note that the outer narrow band is empty when the whole border of the
      // grid is frozen.
      outside_narrow_band_indices.insert(end(outside_narrow_band_indices),
                                         begin(outer_narrow_band_indices),
                                         std::end(outer_narrow_band_indices));

      // Inner dilation bands cannot overlap.
      for (auto k = std::size_t{1}; k < dilation_band_areas.size(); ++k) {
        auto const& inner_dilation_band_indices =
            dilation_bands[dilation_band_areas[k].first];
        assert(!inner_dilation_band_indices.empty());
        auto const inner_narrow_band_indices = NarrowBandDilationBandCells(
            inner_dilation_band_indices, boundary_mask_grid);
        assert(!inner_narrow_band_indices.empty());
        assert(none_of(std::begin(inner_narrow_band_indices),
                       std::end(inner_narrow_band_indices),
                       [=](auto const& distance_grid_index) {
                         return !Inside(distance_grid_index, grid_size);
                       }));
        inside_narrow_band_indices.insert(end(inside_narrow_band_indices),
                                          begin(inner_narrow_band_indices),
                                          std::end(inner_narrow_band_indices));
      }
    }
  }

  return {outside_narrow_band_indices, inside_narrow_band_indices};
}

//! Returns true if the (distance) value @a d  is considered frozen,
//! otherwise false.
//!
//! Preconditions:
//! - @a d is not NaN.
template <typename T>
bool Frozen(T const d) {
  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(!std::isnan(d));
  return -std::numeric_limits<T>::max() < d &&
         d < std::numeric_limits<T>::max();
}

//! Set boundary times on @a time_grid. Times are multiplied by
//! @a multiplier (typically 1 or -1).
//!
//! Preconditions:
//! - Sizes of @a boundary_indices and @a boundary_times are equal.
//!
//! Throws std::invalid_argument if:
//! - The @a check_duplicate_indices is true and there is one or more
//!   duplicate in @a indices.
//! - Not every element in @a boundary_indices is inside @a time_grid.
template <typename T, std::size_t N>
void SetBoundaryCondition(
    std::vector<std::array<std::int32_t, N>> const& boundary_indices,
    std::vector<T> const& boundary_times, T const multiplier,
    bool const check_duplicate_indices, Grid<T, N>* const time_grid) {
  assert(time_grid != nullptr);
  assert(boundary_indices.size() == boundary_times.size() && "Precondition");

  for (auto i = std::size_t{0}; i < boundary_indices.size(); ++i) {
    auto const index = boundary_indices[i];
    auto const time = multiplier * boundary_times[i];
    assert(Inside(index, time_grid->size()) && "Precondition");

    auto& time_cell = time_grid->Cell(index);
    if (check_duplicate_indices) {
      ThrowIfDuplicateBoundaryIndex(Frozen(time_cell), index);
    }
    time_cell = time;
    assert(Frozen(time_cell));
  }
}

//! Returns a (non-null) non-empty narrow band store containing estimated
//! distances for the cells in @a narrow_band_indices. Note that
//! @a narrow_band_indices may contain duplicates.
//!
//! Preconditions:
//! - Boundary condition distances have been set in @a time_grid.
//! - List of narrow band indices is not empty.
//! - Narrow band indices are inside @a time_grid.
//! - Narrow band indices are not frozen in @a time_grid.
template <typename T, std::size_t N, typename E>
std::unique_ptr<NarrowBandStore<T, N>> InitializedNarrowBand(
    std::vector<std::array<std::int32_t, N>> const& narrow_band_indices,
    Grid<T, N> const& time_grid, E const& eikonal_solver) {
  assert(!narrow_band_indices.empty() && "Precondition");

  auto narrow_band =
      std::unique_ptr<NarrowBandStore<T, N>>(new NarrowBandStore<T, N>());
  for (auto const& narrow_band_index : narrow_band_indices) {
    assert(Inside(narrow_band_index, time_grid.size()) && "Precondition");
    assert(!Frozen(time_grid.Cell(narrow_band_index)) && "Precondition");
    narrow_band->Push({eikonal_solver.Solve(narrow_band_index, time_grid),
                       narrow_band_index});
  }
  assert(!narrow_band->empty());

  return narrow_band;
}

//! Compute arrival times using the @a eikonal_solver for the face-neighbors of
//! the cell at @a index. The arrival times are not written to the @a time_grid,
//! but are instead stored in the @a narrow_band.
template <typename T, std::size_t N, typename E>
void UpdateNeighbors(std::array<std::int32_t, N> const& index,
                     E const& eikonal_solver, Grid<T, N>* const time_grid,
                     NarrowBandStore<T, N>* const narrow_band) {
  static_assert(N > 0, "dimensionality cannot be zero");
  static_assert(N == E::kDimension, "mismatching eikonal solver dimension");

  assert(time_grid != nullptr);
  assert(narrow_band != nullptr);
  assert(Inside(index, time_grid->size()));
  assert(Frozen(time_grid->Cell(index)));

  // Update the narrow band. Check face-neighbors in all dimensions.
  auto const kNeighborOffsets = std::array<int32_t, 2>{{-1, 1}};
  for (auto i = std::size_t{0}; i < N; ++i) {
    for (auto const neighbor_offset : kNeighborOffsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;

      if (Inside(neighbor_index, time_grid->size())) {
        // If the neighbor is not frozen compute a distance for it.
        // Note that we don't check if there is an entry for this index
        // in the narrow band already. If we happen to insert multiple
        // distances for the same index the smallest one will be frozen first
        // when marching and the larger distances will be ignored.
        auto& distance_cell = time_grid->Cell(neighbor_index);
        if (!Frozen(distance_cell)) {
          narrow_band->Push({eikonal_solver.Solve(neighbor_index, *time_grid),
                             neighbor_index});
        }
      }
    }
  }
}

//! Compute distances using @a eikonal_solver for all non-frozen cells in
//! @a distance_grid that have a face-connected path to at least one of the
//! cells in @a narrow_band.
//!
//! Preconditions:
//! - @a narrow_band is not empty.
template <typename T, std::size_t N, typename E>
void MarchNarrowBand(E const& eikonal_solver,
                     NarrowBandStore<T, N>* const narrow_band,
                     Grid<T, N>* const time_grid) {
  assert(time_grid != nullptr);
  assert(narrow_band != nullptr);
  assert(!narrow_band->empty() && "Precondition");

  while (!narrow_band->empty()) {
    // Take smallest time from the narrow band and freeze it, i.e.
    // write it to the time grid.
    auto const narrow_band_cell = narrow_band->Pop();
    auto const time = narrow_band_cell.first;
    auto const index = narrow_band_cell.second;

    assert(Inside(index, time_grid->size()));
    auto& time_cell = time_grid->Cell(index);

    // Since we allow multiple values for the same cell index in the
    // narrow band it could happen that this grid cell has already been
    // frozen. In that case just ignore subsequent values from the narrow
    // band for that grid cell and move on.
    if (!Frozen(time_cell)) {
      time_cell = time;
      assert(Frozen(time_cell));

      // Update distances for non-frozen face-neighbors of the newly
      // frozen cell.
      UpdateNeighbors(index, eikonal_solver, time_grid, narrow_band);
    }
  }
}

//! DOCS
//!
//!
//! Throws std::invalid_argument if:
//! - Not the same number of @a indices and @a distances, or
//! - @a indices (and @a distances) are empty, or
//! - Any index is outside the @a distance_grid, or
//! - Any duplicate in @a indices, or
//! - Any value in @a distances does not pass the @a distance_predicate test.
template <typename T, std::size_t N, typename EikonalSolverType, typename P>
std::vector<T> ArrivalTime(
    std::array<std::size_t, N> const& grid_size,
    std::vector<std::array<std::int32_t, N>> const& boundary_indices,
    std::vector<T> const& boundary_times,
    EikonalSolverType const& eikonal_solver, P const boundary_time_predicate,
    bool const negative_inside) {
  typedef T TimeType;

  static_assert(N >= 2, "dimensions must be >= 2");
  static_assert(N == EikonalSolverType::kDimension,
                "mismatching eikonal solver dimension");

  // Check input.
  ThrowIfZeroElementInSize(grid_size);
  ThrowIfEmptyBoundaryIndices(boundary_indices);
  ThrowIfFullGridBoundaryIndices(boundary_indices, grid_size);
  ThrowIfBoundaryIndicesTimesSizeMismatch(boundary_indices, boundary_times);
  std::for_each(std::begin(boundary_indices), std::end(boundary_indices),
                [=](auto const& boundary_index) {
                  ThrowIfBoundaryIndexOutsideGrid(boundary_index, grid_size);
                });
  std::for_each(std::begin(boundary_times), std::end(boundary_times),
                [=](auto const& boundary_time) {
                  ThrowIfInvalidBoundaryTime(
                      boundary_time_predicate(boundary_time), boundary_time);
                });

  auto narrow_band_indices =
      OutsideInsideNarrowBandIndices(boundary_indices, grid_size);
  auto const& outside_narrow_band_indices = narrow_band_indices.first;
  auto const& inside_narrow_band_indices = narrow_band_indices.second;

  auto time_buffer = std::vector<TimeType>(
      LinearSize(grid_size), std::numeric_limits<TimeType>::max());
  assert(std::none_of(std::begin(time_buffer), std::end(time_buffer),
                      [](TimeType const t) { return detail::Frozen(t); }));
  auto time_grid = Grid<TimeType, N>(grid_size, time_buffer);

  if (!inside_narrow_band_indices.empty()) {
    // Set boundaries for marching inside. Always check for duplicate indices.
    auto const check_duplicate_indices = true;
    SetBoundaryCondition(
        boundary_indices, boundary_times,
        TimeType{-1},  // Multiplier, negate boundary times for inside.
        check_duplicate_indices, &time_grid);

    // Initialize inside narrow band with negated boundary times.
    auto inside_narrow_band = InitializedNarrowBand(inside_narrow_band_indices,
                                                    time_grid, eikonal_solver);
    MarchNarrowBand(eikonal_solver, inside_narrow_band.get(), &time_grid);

    if (negative_inside) {
      // Negate all the inside times. Essentially, negate everything
      // computed so far. Note that this also affects the boundary cells.
      std::for_each(std::begin(time_buffer), std::end(time_buffer),
                    [](auto& t) { t = Frozen(t) ? t* TimeType{-1} : t; });
    }
  }

  if (!outside_narrow_band_indices.empty()) {
    // Set boundaries for marching outside. Only check for duplicate indices
    // if this was not done already, i.e. if we marched an inside narrow band.
    auto const check_duplicate_indices = inside_narrow_band_indices.empty();
    SetBoundaryCondition(
        boundary_indices, boundary_times,
        TimeType{1},  // Multiplier, original boundary distances for outside.
        check_duplicate_indices, &time_grid);

    // Initialize outside narrow band with original boundary times.
    auto outside_narrow_band = InitializedNarrowBand(
        outside_narrow_band_indices, time_grid, eikonal_solver);
    MarchNarrowBand(eikonal_solver, outside_narrow_band.get(), &time_grid);
  }

  assert(all_of(std::begin(time_buffer), std::end(time_buffer),
                [](TimeType const t) { return Frozen(t); }));

  return time_buffer;
}

//! Polynomial coefficients are equivalent to std::array index,
//! i.e. Sum(q[i] * x^i) = 0, for i in [0, 2], or simpler
//! q[0] + q[1] * x + q[2] * x^2 = 0.
//!
//! Returns the largest real root, if any (otherwise a NaN value).
//! Note that we are not checking for errors here.
template <typename T>
T SolveEikonalQuadratic(std::array<T, 3> const& q) {
  static_assert(std::is_floating_point<T>::value,
                "quadratic coefficients must be floating point");

  // No error-checking here, caller handles bad values.
  auto const discriminant = q[1] * q[1] - T{4} * q[2] * q[0];
  auto const pos_root = (-q[1] + std::sqrt(discriminant)) / (T{2} * q[2]);
  return pos_root;
}

//! Solve the eikonal equation to get the arrival time (which is distance when
//! @a speed is one) at @a index.
//!
//! The returned value is guaranteed to be positive.
//!
//! Preconditions:
//! - @a speed must be greater than zero.
//! - All elements of @a grid_spacing must be greater than zero.
//! - @a index is inside @a distance_grid.
//! - The cell at @a index must not be frozen in @a distance_grid.
//! - There must be at least one cell in @a distance_grid that is a
//!   frozen face-neighbor of @a index.
//! - Cells in @a distance_grid that are not frozen must have the value
//!   std::numeric_limits<T>::max().
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
T SolveEikonal(std::array<std::int32_t, N> const& index,
               Grid<T, N> const& distance_grid, T const speed,
               std::array<T, N> const& grid_spacing) {
  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(ValidSpeed(speed) && "Precondition");
  assert(ValidGridSpacing(grid_spacing) && "Precondition");
  assert(Inside(index, distance_grid.size()) && "Precondition");
  assert(!Frozen(distance_grid.Cell(index)) && "Precondition");

  // Find the smallest frozen neighbor (if any) in each dimension.
  auto frozen_neighbor_distances = std::array<std::pair<T, std::size_t>, N>();
  auto frozen_neighbor_distances_count = std::size_t{0};
  for (auto i = std::size_t{0}; i < N; ++i) {
    auto neighbor_min_distance = std::numeric_limits<T>::max();
    assert(!Frozen(neighbor_min_distance));

    // Find the smallest face neighbor for this dimension.
    auto neighbor_index = index;

    // -1
    neighbor_index[i] -= int32_t{1};
    if (Inside(neighbor_index, distance_grid.size())) {
      // Note that if the neighbor is not frozen it will have the default
      // distance std::numeric_limits<T>::max().
      auto const neighbor_distance = distance_grid.Cell(neighbor_index);
      if (neighbor_distance < neighbor_min_distance) {
        neighbor_min_distance = neighbor_distance;
        assert(Frozen(neighbor_min_distance));
      }
    }

    // +1
    neighbor_index[i] += int32_t{2};  // -1 + 2 = 1
    if (Inside(neighbor_index, distance_grid.size())) {
      // Note that if the neighbor is not frozen it will have the default
      // distance std::numeric_limits<T>::max().
      auto const neighbor_distance = distance_grid.Cell(neighbor_index);
      if (neighbor_distance < neighbor_min_distance) {
        neighbor_min_distance = neighbor_distance;
        assert(Frozen(neighbor_min_distance));
      }
    }

    // If no frozen neighbor was found that dimension does not contribute
    // to the arrival time.
    if (neighbor_min_distance < std::numeric_limits<T>::max()) {
      frozen_neighbor_distances[frozen_neighbor_distances_count++] = {
          neighbor_min_distance, i};
    }
  }
  assert(frozen_neighbor_distances_count > std::size_t{0} && "Precondition");

  // Define and intiailise q out of if/else for error reporting convenience
  std::array<T, 3> q = std::array<T, 3>{{T{0}, T{0}, T{0}}};

  auto arrival_time = std::numeric_limits<T>::quiet_NaN();
  if (frozen_neighbor_distances_count == 1) {
    // If frozen neighbor in only one dimension we don't need to solve a
    // quadratic.
    auto const distance = frozen_neighbor_distances[0].first;
    auto const j = frozen_neighbor_distances[0].second;
    arrival_time = distance + grid_spacing[j] / speed;
  } else {
    // Initialize quadratic coefficients.
    // auto q = std::array<T, 3>{{T{-1} / Squared(speed), T{0}, T{0}}};
    q[0] = T{-1} / Squared(speed);
    auto const inverse_squared_grid_spacing = InverseSquared(grid_spacing);
    for (auto i = std::size_t{0}; i < frozen_neighbor_distances_count; ++i) {
      auto const distance = frozen_neighbor_distances[i].first;
      auto const j = frozen_neighbor_distances[i].second;
      auto const alpha = inverse_squared_grid_spacing[j];
      q[0] += Squared(distance) * alpha;
      q[1] += T{-2} * distance * alpha;
      q[2] += alpha;
    }
    arrival_time = SolveEikonalQuadratic(q);

    // Fallback when arrival time is negative or NaN.
    if (use_eikonal_fallback && !(arrival_time >= T{0})) {
      // In case the discriminant is negative, we revert to the smallest
      // distance to a single neighbor
      // TODO if dimension N>=3 we could try to find the smallest dimension
      // for the Eikonal equations in dimension N-1
      auto distance = frozen_neighbor_distances[0].first;
      auto j = frozen_neighbor_distances[0].second;
      arrival_time = distance + grid_spacing[j] / speed;
      for (auto i = std::size_t{1}; i < frozen_neighbor_distances_count; ++i) {
        distance = frozen_neighbor_distances[i].first;
        j = frozen_neighbor_distances[i].second;
        arrival_time =
            std::min(arrival_time, distance + grid_spacing[j] / speed);
      }
    }
  }

  ThrowIfInvalidArrivalTimeWithDetails(
      arrival_time, index, speed, q, distance_grid.Cell(index),
      frozen_neighbor_distances, frozen_neighbor_distances_count);
  return arrival_time;
}

//! Solve the eikonal equation to get the arrival time (which is distance when
//! @a speed is one) at @a index. Uses second order derivatives where possible,
//! this version is slower but more accurate. However, it is important to
//! have good boundary conditions that allow second order derivatives to be
//! used early when marching. Otherwise early errors will be propagated.
//!
//! The returned value is guaranteed to be positive.
//!
//! Preconditions:
//! - @a speed must be greater than zero.
//! - All elements of @a grid_spacing must be greater than zero.
//! - @a index is inside @a distance_grid.
//! - The cell at @a index must not be frozen in @a distance_grid.
//! - There must be at least one cell in @a distance_grid that is a
//!   frozen face-neighbor of @a index.
//! - Cells in @a distance_grid that are not frozen must have the value
//!   std::numeric_limits<T>::max().
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
T HighAccuracySolveEikonal(std::array<std::int32_t, N> const& index,
                           Grid<T, N> const& distance_grid, T const speed,
                           std::array<T, N> const& grid_spacing) {
  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(ValidSpeed(speed) && "Precondition");
  assert(ValidGridSpacing(grid_spacing) && "Precondition");
  assert(Inside(index, distance_grid.size()) && "Precondition");
  assert(!Frozen(distance_grid.Cell(index)) && "Precondition");

  // Find the smallest frozen neighbor(s) (if any) in each dimension.
  auto const neighbor_offsets = std::array<int32_t, 2>{{-1, 1}};
  auto frozen_neighbor_distances =
      std::array<std::pair<std::pair<T, T>, std::size_t>, N>();
  auto frozen_neighbor_distances_count = std::size_t{0};
  for (auto i = std::size_t{0}; i < N; ++i) {
    auto neighbor_min_distance = std::numeric_limits<T>::max();
    auto neighbor_min_distance2 = std::numeric_limits<T>::max();
    assert(!Frozen(neighbor_min_distance));
    assert(!Frozen(neighbor_min_distance2));

    // Check neighbors in both directions for this dimenion.
    for (auto const neighbor_offset : neighbor_offsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;
      if (Inside(neighbor_index, distance_grid.size())) {
        auto const neighbor_distance = distance_grid.Cell(neighbor_index);
        if (neighbor_distance < neighbor_min_distance) {
          // Neighbor one step away is frozen.
          assert(Frozen(neighbor_distance));
          neighbor_min_distance = neighbor_distance;

          // Check if neighbor two steps away is frozen and has smaller
          // (or equal) distance than neighbor one step away. Reset
          // the distance first since otherwise we might get the secondary
          // distance from the previous neighbor offset.
          neighbor_min_distance2 = std::numeric_limits<T>::max();
          auto neighbor_index2 = neighbor_index;
          neighbor_index2[i] += neighbor_offset;
          if (Inside(neighbor_index2, distance_grid.size())) {
            auto const neighbor_distance2 = distance_grid.Cell(neighbor_index2);
            if (neighbor_distance2 <= neighbor_distance) {
              // Neighbor index two steps away is frozen.
              assert(Frozen(neighbor_distance2));
              neighbor_min_distance2 = neighbor_distance2;
            }
          }
        }
      }
    }

    if (neighbor_min_distance2 < std::numeric_limits<T>::max()) {
      // Two frozen neighbors in this dimension.
      assert(neighbor_min_distance < std::numeric_limits<T>::max());
      frozen_neighbor_distances[frozen_neighbor_distances_count++] = {
          {neighbor_min_distance, neighbor_min_distance2}, i};
    } else if (neighbor_min_distance < std::numeric_limits<T>::max()) {
      // One frozen neighbor in this dimension.
      frozen_neighbor_distances[frozen_neighbor_distances_count++] = {
          {neighbor_min_distance, std::numeric_limits<T>::max()}, i};
    }
    // else: no frozen neighbors in this dimension.
  }
  assert(frozen_neighbor_distances_count > std::size_t{0} && "Precondition");

  // Define and intiailise q out of if/else for error reporting convenience
  std::array<T, 3> q = std::array<T, 3>{{T{0}, T{0}, T{0}}};

  auto arrival_time = std::numeric_limits<T>::quiet_NaN();
  if (frozen_neighbor_distances_count == 1) {
    // If frozen neighbor in only one dimension we don't need to solve a
    // quadratic.
    auto const distance = frozen_neighbor_distances[0].first.first;
    auto const j = frozen_neighbor_distances[0].second;
    arrival_time = distance + grid_spacing[j] / speed;
  } else {
    // Initialize quadratic coefficients.
    // auto q = std::array<T, 3>{{T{-1} / Squared(speed), T{0}, T{0}}};
    q[0] = T{-1} / Squared(speed);
    auto const inverse_squared_grid_spacing = InverseSquared(grid_spacing);

    for (auto i = std::size_t{0}; i < frozen_neighbor_distances_count; ++i) {
      auto const distance = frozen_neighbor_distances[i].first.first;
      auto const distance2 = frozen_neighbor_distances[i].first.second;
      auto const j = frozen_neighbor_distances[i].second;
      if (distance2 < std::numeric_limits<T>::max()) {
        // Second order coefficients.
        assert(distance < std::numeric_limits<T>::max());
        auto const alpha = (T{9} / T{4}) * inverse_squared_grid_spacing[j];
        auto const t = (T{1} / T{3}) * (T{4} * distance - distance2);
        q[0] += Squared(t) * alpha;
        q[1] += T{-2} * t * alpha;
        q[2] += alpha;
      } else if (distance < std::numeric_limits<T>::max()) {
        // First order coefficients.
        auto const alpha = inverse_squared_grid_spacing[j];
        q[0] += Squared(distance) * alpha;
        q[1] += T{-2} * distance * alpha;
        q[2] += alpha;
      }
    }
    arrival_time = SolveEikonalQuadratic(q);

    // Fallback when arrival time is negative or NaN.
    if (use_eikonal_fallback && !(arrival_time >= T{0})) {
      // In case the discriminant is negative, we revert to the smallest
      // distance to a single neighbor
      // TODO if dimension N>=3 we could try to find the smallest dimension
      // for the Eikonal equations in dimension N-1
      auto distance = frozen_neighbor_distances[0].first.first;
      auto j = frozen_neighbor_distances[0].second;
      arrival_time = distance + grid_spacing[j] / speed;
      for (auto i = std::size_t{1}; i < frozen_neighbor_distances_count; ++i) {
        distance = frozen_neighbor_distances[i].first.first;
        j = frozen_neighbor_distances[i].second;
        arrival_time =
            std::min(arrival_time, distance + grid_spacing[j] / speed);
      }
    }
  }

  ThrowIfInvalidArrivalTimeWithHighAccuracyDetails(
      arrival_time, index, speed, q, distance_grid.Cell(index),
      frozen_neighbor_distances, frozen_neighbor_distances_count);
  return arrival_time;
}

//!
//!
//! Implementation follows pseudo-code given in "Fluid Simulation for
//! Computer Graphics" by Robert Bridson.
//!
//! Preconditions:
//! - @a dx must be greater than zero.
//! - @a index is inside @a distance_grid.
//! - The cell at @a index must not be frozen in @a distance_grid.
//! - There must be at least one cell in @a distance_grid that is a
//!   frozen face-neighbor of @a index.
//! - Cells in @a distance_grid that are not frozen must have the value
//!   std::numeric_limits<T>::max().
//!
//! Note: Currently supports only uniformly spaced (square) cells.
//! Note: Currently supports only 1D, 2D, and 3D.
template <typename T, std::size_t N>
T SolveDistance(std::array<std::int32_t, N> const& index,
                Grid<T, N> const& distance_grid, T const dx) {
  static_assert(1 <= N && N <= 3, "invalid dimensionality");

  assert(Inside(index, distance_grid.size()) && "Precondition");
  assert(!Frozen(distance_grid.Cell(index)) && "Precondition");

  auto phi = std::array<T, N>();
  std::fill(std::begin(phi), std::end(phi), std::numeric_limits<T>::max());
  auto phi_count = std::size_t{0};

  // Find the smallest frozen neighbor(s) (if any) in each dimension.
  for (auto i = std::size_t{0}; i < N; ++i) {
    auto neighbor_min_distance = std::numeric_limits<T>::max();
    assert(!Frozen(neighbor_min_distance));

    // -1
    auto neighbor_index = index;
    neighbor_index[i] -= 1;
    if (Inside(neighbor_index, distance_grid.size())) {
      auto const neighbor_distance = distance_grid.Cell(neighbor_index);
      if (neighbor_distance < neighbor_min_distance) {
        neighbor_min_distance = neighbor_distance;
      }
    }

    // -1 + 2 = +1
    neighbor_index[i] += 2;
    if (Inside(neighbor_index, distance_grid.size())) {
      auto const neighbor_distance = distance_grid.Cell(neighbor_index);
      if (neighbor_distance < neighbor_min_distance) {
        neighbor_min_distance = neighbor_distance;
      }
    }

    if (neighbor_min_distance < std::numeric_limits<T>::max()) {
      phi[phi_count++] = neighbor_min_distance;
    }
  }
  assert(phi_count > 0 && "Precondition");

  // Sort ascending using a sorting network approach.
  if (N >= 2 && phi[0] > phi[1]) {
    std::swap(phi[0], phi[1]);
  }
  if (N == 3 && phi[1] > phi[2]) {
    std::swap(phi[1], phi[2]);
  }
  if (N == 3 && phi[0] > phi[1]) {
    std::swap(phi[0], phi[1]);
  }

  auto distance = phi[0] + dx;
  if (N >= 2 && phi_count > 1 && distance > phi[1]) {
    distance =
        T(0.5) * (phi[0] + phi[1] +
                  std::sqrt(T(2) * Squared(dx) - Squared(phi[1] - phi[0])));
    if (N == 3 && phi_count == 3 && distance > phi[2]) {
      auto const phi_sum = phi[0] + phi[1] + phi[2];
      auto phi_sum_squared =
          Squared(phi[0]) + Squared(phi[1]) + Squared(phi[2]);
      distance =
          (T(1) / T(3)) *
          (phi_sum +
           std::sqrt(std::max(T(0), Squared(phi_sum) - T(3) * (phi_sum_squared -
                                                               Squared(dx)))));
    }
  }

  // Arrival time is distance here since we have assumed that speed is one.
  ThrowIfInvalidArrivalTime(distance, index);
  return distance;
}

//! Base class for Eikonal solvers.
//! Note: dtor is not virtual!
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
class EikonalSolverBase {
 public:
  typedef T ScalarType;
  static std::size_t const kDimension = N;
  static bool const allowEikonalFallback = use_eikonal_fallback;

 protected:
  explicit EikonalSolverBase(std::array<T, N> const& grid_spacing)
      : grid_spacing_(grid_spacing) {
    ThrowIfInvalidGridSpacing(grid_spacing_);
  }

  std::array<T, N> const& grid_spacing() const { return grid_spacing_; }

 private:
  std::array<T, N> const grid_spacing_;
};

//! Base class for Eikonal solvers with uniform speed.
//! Note: dtor is not virtual!
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
class UniformSpeedEikonalSolverBase
    : public EikonalSolverBase<T, N, use_eikonal_fallback> {
 protected:
  UniformSpeedEikonalSolverBase(std::array<T, N> const& grid_spacing,
                                T const uniform_speed)
      : EikonalSolverBase<T, N, use_eikonal_fallback>(grid_spacing),
        uniform_speed_(uniform_speed) {
    ThrowIfZeroOrNegativeOrNanSpeed(uniform_speed_);
  }

  //! Returns the uniform speed, guaranteed to be:
  //! - Non-zero
  //! - Positive
  //! - Not NaN
  T uniform_speed() const { return uniform_speed_; }

 private:
  T const uniform_speed_;
};

//! Base class for Eikonal solvers with varying speed.
//! Note: dtor is not virtual!
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
class VaryingSpeedEikonalSolverBase
    : public EikonalSolverBase<T, N, use_eikonal_fallback> {
 protected:
  VaryingSpeedEikonalSolverBase(
      std::array<T, N> const& grid_spacing,
      std::array<std::size_t, N> const& speed_grid_size,
      std::vector<T> const& speed_buffer)
      : EikonalSolverBase<T, N, use_eikonal_fallback>(grid_spacing),
        speed_grid_(speed_grid_size, speed_buffer) {
    for (auto const speed : speed_buffer) {
      ThrowIfZeroOrNegativeOrNanSpeed(speed);
    }
  }

  //! Returns the speed at @a index in the speed grid, guaranteed to be:
  //! - Non-zero
  //! - Positive
  //! - Not NaN
  //!
  //! Throws an std::invalid_argument exception if @a index is outside the
  //! speed grid.
  T Speed(std::array<std::int32_t, N> const& index) const {
    ThrowIfSpeedIndexOutsideGrid(index, speed_grid_.size());
    return speed_grid_.Cell(index);
  }

 private:
  ConstGrid<T, N> const speed_grid_;
};

}  // namespace detail

//! Provides methods for solving the eikonal equation for a single grid cell
//! at a time using the current distance grid. Uses a uniform speed for
//! the entire grid.
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
class UniformSpeedEikonalSolver
    : public detail::UniformSpeedEikonalSolverBase<T, N, use_eikonal_fallback> {
 public:
  explicit UniformSpeedEikonalSolver(std::array<T, N> const& grid_spacing,
                                     T const uniform_speed = T{1})
      : detail::UniformSpeedEikonalSolverBase<T, N, use_eikonal_fallback>(
            grid_spacing, uniform_speed) {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) of other cells.
  T Solve(std::array<std::int32_t, N> const& index,
          detail::Grid<T, N> const& distance_grid) const {
    return detail::SolveEikonal<T, N, use_eikonal_fallback>(
        index, distance_grid, this->uniform_speed(), this->grid_spacing());
  }
};

//! Provides methods for solving the eikonal equation for a single grid cell
//! at a time using the current distance grid. Uses a uniform speed for
//! the entire grid. When possible uses second order derivates to achieve
//! better accuracy.
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
class HighAccuracyUniformSpeedEikonalSolver
    : public detail::UniformSpeedEikonalSolverBase<T, N, use_eikonal_fallback> {
 public:
  explicit HighAccuracyUniformSpeedEikonalSolver(
      std::array<T, N> const& grid_spacing, T const uniform_speed = T{1})
      : detail::UniformSpeedEikonalSolverBase<T, N, use_eikonal_fallback>(
            grid_spacing, uniform_speed) {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) of other cells.
  T Solve(std::array<std::int32_t, N> const& index,
          detail::Grid<T, N> const& distance_grid) const {
    return detail::HighAccuracySolveEikonal<T, N, use_eikonal_fallback>(
        index, distance_grid, this->uniform_speed(), this->grid_spacing());
  }
};

//! Provides methods for solving the eikonal equation for a single grid cell
//! at a time using the current distance grid. A speed grid must be provided
//! and that grid must cover the arrival time grid.
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
class VaryingSpeedEikonalSolver
    : public detail::VaryingSpeedEikonalSolverBase<T, N, use_eikonal_fallback> {
 public:
  VaryingSpeedEikonalSolver(std::array<T, N> const& grid_spacing,
                            std::array<std::size_t, N> const& speed_grid_size,
                            std::vector<T> const& speed_buffer)
      : detail::VaryingSpeedEikonalSolverBase<T, N, use_eikonal_fallback>(
            grid_spacing, speed_grid_size, speed_buffer) {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) of other cells.
  T Solve(std::array<std::int32_t, N> const& index,
          detail::Grid<T, N> const& distance_grid) const {
    return detail::SolveEikonal<T, N, use_eikonal_fallback>(
        index, distance_grid, this->Speed(index), this->grid_spacing());
  }
};

//! Provides methods for solving the eikonal equation for a single grid cell
//! at a time using the current distance grid. A speed grid must be provided
//! and that grid must cover the arrival time grid. When possible uses second
//! order derivates to achieve better accuracy.
template <typename T, std::size_t N, bool use_eikonal_fallback = true>
class HighAccuracyVaryingSpeedEikonalSolver
    : public detail::VaryingSpeedEikonalSolverBase<T, N, use_eikonal_fallback> {
 public:
  HighAccuracyVaryingSpeedEikonalSolver(
      std::array<T, N> const& grid_spacing,
      std::array<std::size_t, N> const& speed_grid_size,
      std::vector<T> const& speed_buffer)
      : detail::VaryingSpeedEikonalSolverBase<T, N, use_eikonal_fallback>(
            grid_spacing, speed_grid_size, speed_buffer) {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) of other cells.
  T Solve(std::array<std::int32_t, N> const& index,
          detail::Grid<T, N> const& distance_grid) const {
    return detail::HighAccuracySolveEikonal<T, N, use_eikonal_fallback>(
        index, distance_grid, this->Speed(index), this->grid_spacing());
  }
};

//! Provides methods for solving the eikonal equation for a single grid cell
//! at a time using the current distance grid. The speed is assumed to be
//! one for the entire grid, meaning that arrival time can be interpreted
//! as distance.
//!
//! Note: Currently only supports uniform grid spacing.
template <typename T, std::size_t N>
class DistanceSolver {
 public:
  typedef T ScalarType;
  static std::size_t const kDimension = N;

  explicit DistanceSolver(T const dx) : dx_(dx) {
    detail::ThrowIfInvalidGridSpacing(detail::FilledArray<T, N>(dx_));
  }

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) of other cells.
  T Solve(std::array<std::int32_t, N> const& index,
          detail::Grid<T, N> const& distance_grid) const {
    return detail::SolveDistance(index, distance_grid, dx_);
  }

 private:
  T const dx_;
};

//! Compute the signed distance on a grid.
//!
//! Input:
//!   grid_size          - Number of grid cells in each dimension.
//!   boundary_indices   - Integer coordinates of cells with provided distances.
//!   boundary_distances - Signed distances assigned to boundary cells.
//!
//! Preconditions:
//!   - grid_size may not have a zero element.
//!   - frozen_indices, frozen_distances and normals must have the same size.
//!   - frozen_indices must all be within size.
//!
//! TODO - example usage!
template <typename T, std::size_t N, typename EikonalSolverType>
std::vector<T> SignedArrivalTime(
    std::array<std::size_t, N> const& grid_size,
    std::vector<std::array<std::int32_t, N>> const& boundary_indices,
    std::vector<T> const& boundary_times,
    EikonalSolverType const& eikonal_solver) {
  auto const boundary_time_predicate = [](auto const t) {
    return !std::isnan(t) && detail::Frozen(t);
  };
  auto constexpr negative_inside = true;
  return detail::ArrivalTime(grid_size, boundary_indices, boundary_times,
                             eikonal_solver, boundary_time_predicate,
                             negative_inside);
}

}  // namespace fast_marching_method
}  // namespace thinks

#endif  // INCLUDE_THINKS_FAST_MARCHING_METHOD_FAST_MARCHING_METHOD_HPP_
