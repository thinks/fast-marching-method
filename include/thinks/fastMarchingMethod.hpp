#ifndef THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED
#define THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <stack>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>


namespace thinks {
namespace fast_marching_method {
namespace detail {

//! Returns the product of the elements in array @a size.
//! Throws a std::invalid_argument exception if one or more of the elements
//! in @a size are zero.
//!
//! Note: Not checking for integer overflow here!
template<std::size_t N> inline
std::size_t LinearSize(std::array<std::size_t, N> const& size)
{
  using namespace std;

  ThrowIfZeroElementInSize(size);
  return accumulate(begin(size), end(size), size_t{1}, multiplies<size_t>());
}


//! Returns x * x.
template<typename T> inline constexpr
T Squared(T const x)
{
  return x * x;
}


//! Returns 1 / (x * x).
template<typename T> inline constexpr
T InverseSquared(T const x)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "scalar type must be floating point");

  return T{1} / Squared(x);
}

//! Returns element-wise 1 / (a * a).
template<typename T, std::size_t N> inline
std::array<T, N> InverseSquared(std::array<T, N> const& a)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "scalar type must be floating point");

  auto r = array<T, N>();
  transform(begin(a), end(a), begin(r),
            [](T const x) { return InverseSquared(x); });
  return r;
}


//! Returns true if @a index is inside @a size, otherwise false.
template<std::size_t N> inline
bool Inside(
  std::array<std::int32_t, N> const& index,
  std::array<std::size_t, N> const& size)
{
  using namespace std;

  for (auto i = size_t{0}; i < N; ++i) {
    // Cast is safe since we check that index[i] is greater than or
    // equal to zero first.
    if (!(0 <= index[i] && static_cast<size_t>(index[i]) < size[i])) {
      return false;
    }
  }

  return true;
}


//! Returns a string representation of the array @a a.
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


//! Throws a std::invalid_argument if one or more of the elements in @a size
//! is zero.
template<std::size_t N> inline
void ThrowIfZeroElementInSize(std::array<std::size_t, N> const& size)
{
  using namespace std;

  if (find_if(begin(size), end(size),
              [](auto const x) { return x == size_t{0}; }) != end(size)) {
    auto ss = stringstream();
    ss << "invalid size: " << ToString(size);
    throw invalid_argument(ss.str());
  }
}


//! Throws a std::invalid_argument if the linear size of @a grid_size is not
//! equal to @a cell_buffer_size.
template<std::size_t N> inline
void ThrowIfInvalidCellBufferSize(
  std::array<std::size_t, N> const& grid_size,
  std::size_t const cell_buffer_size)
{
  using namespace std;

  if (LinearSize(grid_size) != cell_buffer_size) {
    auto ss = stringstream();
    ss << "grid size " << ToString(grid_size)
       << " does not match cell buffer size " << cell_buffer_size;
    throw invalid_argument(ss.str());
  }
}


//! Throws a std::invalid_argument if one or more of the elements in
//! @a grid_spacing is less than or equal to zero (or NaN).
template<typename T, std::size_t N> inline
void ThrowIfInvalidGridSpacing(std::array<T, N> const& grid_spacing)
{
  using namespace std;

  if (find_if(begin(grid_spacing), end(grid_spacing),
              [](auto const x) { return isnan(x) || x <= T(0); }) !=
      end(grid_spacing)) {
    auto ss = stringstream();
    ss << "invalid grid spacing: " << ToString(grid_spacing);
    throw invalid_argument(ss.str());
  }
}


//! Throws a std::invalid_argument if @a speed is less than or equal
//! to zero or NaN.
template<typename T> inline
void ThrowIfZeroOrNegativeOrNanSpeed(T const speed)
{
  using namespace std;

  if (isnan(speed) || speed <= T{0}) {
    auto ss = stringstream();
    ss << "invalid speed: " << speed;
    throw invalid_argument(ss.str());
  }
}


//! Throw a std::invalid_argument if @a boundary_indices is empty.
template<std::size_t N> inline
void ThrowIfEmptyBoundaryIndices(
  std::vector<std::array<std::int32_t, N>> const& boundary_indices)
{
  using namespace std;

  if (boundary_indices.empty()) {
    throw invalid_argument("empty boundary condition");
  }
}


//!
//!
template<std::size_t N> inline
void ThrowIfBoundaryIndexOutsideGrid(
  std::array<std::int32_t, N> const& boundary_index,
  std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  if (!Inside(boundary_index, grid_size)) {
    auto ss = stringstream();
    ss << "boundary index outside grid - "
       << "index: " << ToString(boundary_index) << ", "
       << "grid size: " << ToString(grid_size);
    throw invalid_argument(ss.str());
  }
}


//!
//!
template<std::size_t N> inline
void ThrowIfBoundaryIndexOutsideGrid(
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  for (auto const& boundary_index : boundary_indices) {
    ThrowIfBoundaryIndexOutsideGrid(boundary_index, grid_size);
  }
}


//!
template<typename T> inline
void ThrowIfInvalidBoundaryDistance(bool const valid, T const distance)
{
  using namespace std;

  if (!valid) {
    auto ss = stringstream();
    ss << "invalid boundary distance: " << distance;
    throw invalid_argument(ss.str());
  }
}


//!
template<std::size_t N> inline
void ThrowIfDuplicateBoundaryIndex(
  bool const duplicate,
  std::array<std::int32_t, N> const& index)
{
  using namespace std;

  if (duplicate) {
    auto ss = stringstream();
    ss << "duplicate boundary index: " << ToString(index);
    throw invalid_argument(ss.str());
  }
}


//! Returns an array that can be used to transform an N-dimensional index
//! into a linear index.
template<std::size_t N>
std::array<std::size_t, N - 1> GridStrides(
  std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  auto strides = array<size_t, N - 1>();
  auto stride = size_t{1};
  for (auto i = size_t{1}; i < N; ++i) {
    stride *= grid_size[i - 1];
    strides[i - 1] = stride;
  }
  return strides;
}


//! Returns a linear (scalar) index into an array representing an
//! N-dimensional grid for integer coordinate @a index.
//! Note that this function does not check for integer overflow!
template<std::size_t N> inline
std::size_t GridLinearIndex(
  std::array<std::int32_t, N> const& index,
  std::array<std::size_t, N - 1> const& grid_strides)
{
  using namespace std;

  auto k = static_cast<size_t>(index[0]);
  for (auto i = size_t{1}; i < N; ++i) {
    k += index[i] * grid_strides[i - 1];
  }
  return k;
}


//! Access a linear array as if it were an N-dimensional grid.
//! Allows mutating operations on the underlying array.
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
//!   std::cout << "cell (0,0): " << grid.cell({{0, 0}}) << std::endl;
//!   std::cout << "cell (1,0): " << grid.cell({{1, 0}}) << std::endl;
//!   std::cout << "cell (0,1): " << grid.cell({{0, 1}}) << std::endl;
//!   std::cout << "cell (1,1): " << grid.cell({{1, 1}}) << std::endl;
//!
//!   grid.Cell({{0, 1}}) = 5.3;
//!   std::cout << "cell (0,1): " << grid.cell({{0, 1}}) << std::endl;
//!
//! Output:
//!   cell (0,0): 0
//!   cell (1,0): 1
//!   cell (0,1): 2
//!   cell (1,1): 3
//!
//!   cell (0,1): 5.3
template<typename T, std::size_t N>
class Grid
{
public:
  typedef T CellType;
  typedef std::array<std::size_t, N> SizeType;
  typedef std::array<std::int32_t, N> IndexType;

  //! Construct a grid from a given @a size and @a cell_buffer. Does not
  //! take ownership of the cell buffer, it is assumed that this buffer exists
  //! during the life-time of the grid object.
  //!
  //! Throws a std::invalid_argument if:
  //! - size evaluates to a zero linear size, i.e. if any of the
  //!   elements are zero.
  //! - the linear size of @a size is not equal to the @a cell_buffer size.
  Grid(SizeType const& size, std::vector<T>& cell_buffer)
    : size_(size)
    , strides_(GridStrides(size))
    , cells_(nullptr)
  {
    ThrowIfZeroElementInSize(size);
    ThrowIfInvalidCellBufferSize(size, cell_buffer.size());

    assert(!cell_buffer.empty());
    cells_ = &cell_buffer.front();
  }

  SizeType size() const
  {
    return size_;
  }

  //! Returns a reference to the cell at @a index. No range checking!
  CellType& Cell(IndexType const& index)
  {
    assert(GridLinearIndex(index, strides_) < LinearSize(size()));
    return cells_[GridLinearIndex(index, strides_)];
  }

  //! Returns a const reference to the cell at @a index. No range checking!
  CellType const& Cell(IndexType const& index) const
  {
    assert(GridLinearIndex(index, strides_) < LinearSize(size()));
    return cells_[GridLinearIndex(index, strides_)];
  }

private:
  std::array<std::size_t, N> const size_;
  std::array<std::size_t, N - 1> const strides_;
  CellType* cells_;
};


//! Access a linear array as if it were an N-dimensional grid.
//! Does not allow mutating operations on the underlying array.
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
template<typename T, std::size_t N>
class ConstGrid
{
public:
  typedef T CellType;
  typedef std::array<std::size_t, N> SizeType;
  typedef std::array<std::int32_t, N> IndexType;

  //! Construct a grid from a given @a size and @a cell_buffer. Does not
  //! take ownership of the cell buffer, it is assumed that this buffer exists
  //! during the life-time of the grid object.
  //!
  //! Throws a std::invalid_argument if:
  //! - size evaluates to a zero linear size, i.e. if any of the
  //!   elements are zero.
  //! - the linear size of @a size is not equal to the @a cell_buffer size.
  ConstGrid(SizeType const& size, std::vector<T> const& cell_buffer)
    : size_(size)
    , strides_(GridStrides(size))
    , cells_(nullptr)
  {
    ThrowIfZeroElementInSize(size);
    ThrowIfInvalidCellBufferSize(size, cell_buffer.size());

    assert(!cell_buffer.empty());
    cells_ = &cell_buffer.front();
  }

  SizeType size() const
  {
    return size_;
  }

  //! Returns a const reference to the cell at @a index. No range checking!
  CellType const& Cell(IndexType const& index) const
  {
    assert(GridLinearIndex(index, strides_) < LinearSize(size()));
    return cells_[GridLinearIndex(index, strides_)];
  }

private:
  std::array<std::size_t, N> const size_;
  std::array<std::size_t, N - 1> const strides_;
  CellType const* cells_;
};


//! Acceleration structure for keeping track of the smallest distance cell
//! in the narrow band.
template<typename T, std::size_t N>
class NarrowBandStore
{
public:
  typedef T DistanceType;
  typedef std::array<std::int32_t, N> IndexType;
  typedef std::pair<DistanceType, IndexType> ValueType;

  NarrowBandStore()
  {}

  bool empty() const
  {
    return min_heap_.empty();
  }

  ValueType Pop()
  {
    assert(!min_heap_.empty());
    auto const v = min_heap_.top(); // O(1)
    min_heap_.pop(); // O(log N)
    return v;
  }

  void Push(ValueType const& value)
  {
    min_heap_.push(value); // O(log N)
  }

private:
  struct HeapComparison_
  {
    bool operator()(ValueType const& lhs, ValueType const& rhs)
    {
      return lhs.first > rhs.first;
    }
  };

  typedef std::priority_queue<
    ValueType,
    std::vector<ValueType>,
    HeapComparison_> MinHeap_;

  MinHeap_ min_heap_;
};


template<std::size_t N>
class IndexIterator
{
public:
  explicit IndexIterator(std::array<std::size_t, N> const& size)
    : size_(size)
  {
    using namespace std;

    if (LinearSize(size) == size_t{0}) {
      throw runtime_error("zero size element");
    }

    fill(begin(index_), end(index_), int32_t{0});
  }

  std::array<std::int32_t, N> index() const
  {
    return index_;
  }

  bool Next()
  {
    using namespace std;

    auto i = int32_t{N - 1};
    while (i >= 0) {
      assert(size_[i] > size_t{0});
      if (static_cast<size_t>(index_[i]) < size_[i] - 1) {
        ++index_[i];
        return true;
      }
      else {
        index_[i--] = 0;
      }
    }
    return false;
  }

private:
  std::array<std::size_t, N> const size_;
  std::array<std::int32_t, N> index_;
};


//! Returns base^exponent as a compile-time constant.
//! Note: Not checking for integer overflow here!
constexpr std::size_t inline
static_pow(std::size_t const base, std::size_t const exponent)
{
  using namespace std;

  // NOTE: Cannot use loops in constexpr functions in C++11, have to use
  // recursion here.
  return exponent == size_t{0} ?
    size_t{1} :
    base * static_pow(base, exponent - 1);
}


//! Returns an array of pairs, where each element is the min/max index
//! coordinates in the corresponding dimension.
//!
//! Preconditions:
//! - @a indices is not empty.
template<std::size_t N> inline
std::array<std::pair<std::int32_t, std::int32_t>, N> BoundingBox(
  std::vector<std::array<std::int32_t, N>> const& indices)
{
  using namespace std;

  static_assert(N > 0, "Dimensionality cannot be zero");

  assert(!indices.empty() && "Precondition");

  // Initialize bounding box in all dimensions.
  auto bbox = array<pair<int32_t, int32_t>, N>();
  for (auto i = size_t{0}; i < N; ++i) {
    bbox[i].first = numeric_limits<int32_t>::max();
    bbox[i].second = numeric_limits<int32_t>::min();
  }

  // Update with each index in all dimensions.
  for (auto const& index : indices) {
    for (auto i = size_t{0}; i < N; ++i) {
      bbox[i].first = min(bbox[i].first, index[i]);
      bbox[i].second = max(bbox[i].second, index[i]);
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
template<std::size_t N> inline
std::size_t HyperVolume(
  std::array<std::pair<std::int32_t, std::int32_t>, N> const& bbox)
{
  using namespace std;

  auto hyper_volume = size_t{1};
  for (auto i = size_t{0}; i < N; ++i) {
    assert(bbox[i].first <= bbox[i].second);
    hyper_volume *= (bbox[i].second - bbox[i].first + 1);
  }
  return hyper_volume;
}


//! Returns true if the (distance) value @a d  is considered frozen,
//! otherwise false.
//!
//! Preconditions:
//! - @a d is not NaN.
template <typename T> inline
bool Frozen(T const d)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(!isnan(d));
  return -numeric_limits<T>::max() < d && d < numeric_limits<T>::max();
}


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
  auto valid_index = true;
  while (valid_index) {
    auto offset = index_iter.index();
    for_each(begin(offset), end(offset), [](auto& d) { d -= int32_t{1}; });
    if (!all_of(begin(offset), end(offset), [](auto const i){ return i == 0; })) {
      neighbor_offsets[offset_index++] = offset;
    }
    valid_index = index_iter.Next();
  }
  assert(offset_index == static_pow(3, N) - 1);

  return neighbor_offsets;
}


template<std::size_t N> inline
std::array<std::array<std::int32_t, N>, 2 * N> FaceNeighborOffsets()
{
  using namespace std;

  auto neighbor_offsets = array<array<int32_t, N>, size_t{2} * N>{};
  for (auto i = size_t{0}; i < N; ++i) {
    for (auto j = size_t{0}; j < N; ++j) {
      if (j == i) {
        neighbor_offsets[2 * i + 0][j] = int32_t{+1};
        neighbor_offsets[2 * i + 1][j] = int32_t{-1};
      }
      else {
        neighbor_offsets[2 * i + 0][j] = int32_t{0};
        neighbor_offsets[2 * i + 1][j] = int32_t{0};
      }
    }
  }
  return neighbor_offsets;
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
template<std::size_t N, typename NeighborOffsetIt> inline
std::vector<std::vector<std::array<std::int32_t, N>>>
ConnectedComponents(
  std::vector<std::array<std::int32_t, N>> const& foreground_indices,
  std::array<std::size_t, N> const& grid_size,
  NeighborOffsetIt const neighbor_offset_begin,
  NeighborOffsetIt const neighbor_offset_end)
{
  // TODO: This is a very naive implementation of finding connected
  //       components. It is nice that it generalizes easily to arbitrary
  //       dimensions, but it is very inefficient. Should use the method
  //       that builds equivalence trees which is O(n).

  using namespace std;

  enum class LabelCell : uint8_t {
    kBackground = uint32_t{0},
    kForeground,
    kLabelled
  };

  if (foreground_indices.empty()) {
    return vector<vector<array<int32_t, N>>>();
  }

  auto label_buffer =
    vector<LabelCell>(LinearSize(grid_size), LabelCell::kBackground);
  auto label_grid = Grid<LabelCell, N>(grid_size, label_buffer);

  for (auto const& foreground_index : foreground_indices) {
    // Note: We don't check for duplicate indices here since it doesn't
    //       affect the algorithm.
    assert(Inside(foreground_index, label_grid.size()));
    label_grid.Cell(foreground_index) = LabelCell::kForeground;
  }

  auto connected_components = vector<vector<array<int32_t, N>>>();
  for (auto const& foreground_index : foreground_indices) {
    assert(Inside(foreground_index, label_grid.size()));
    assert(label_grid.Cell(foreground_index) == LabelCell::kForeground ||
           label_grid.Cell(foreground_index) == LabelCell::kLabelled);

    if (label_grid.Cell(foreground_index) == LabelCell::kForeground) {
      // This index has not already been labelled.
      // Start a new component.
      auto component = vector<array<int32_t, N>>();
      auto neighbor_indices = stack<array<int32_t, N>>();
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
          for (auto i = size_t{0}; i < N; ++i) {
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
template<std::size_t N> inline
std::array<std::int32_t, N>
DistanceGridIndexFromDilationGridIndex(
  std::array<std::int32_t, N> const& dilation_grid_index)
{
  auto distance_grid_index = dilation_grid_index;
  for_each(begin(distance_grid_index), end(distance_grid_index),
           [](auto& d) { d -= int32_t{1}; });
  return distance_grid_index;
}


//! Returns @a distance_grid_index transformed to the dilation grid.
template<std::size_t N> inline
std::array<std::int32_t, N>
DilationGridIndexFromDistanceGridIndex(
  std::array<std::int32_t, N> const& distance_grid_index)
{
  auto dilation_grid_index = distance_grid_index;
  for_each(begin(dilation_grid_index), end(dilation_grid_index),
           [](auto& d) { d += int32_t{1}; assert(d >= int32_t{0}); });
  return dilation_grid_index;
}


//!
//!
//! Preconditions:
//! - All elements in @a grid_indices are inside @a grid_size.
//! - Neighbor offsets are not larger than one cell in any dimension.
template<
  std::size_t N,
  typename DilationNeighborOffsetIt,
  typename BandNeighborOffsetIt> inline
std::vector<std::vector<std::array<std::int32_t, N>>>
DilationBands(
  std::vector<std::array<std::int32_t, N>> const& grid_indices,
  std::array<std::size_t, N> const& grid_size,
  DilationNeighborOffsetIt const dilation_neighbor_offset_begin,
  DilationNeighborOffsetIt const dilation_neighbor_offset_end,
  BandNeighborOffsetIt const band_neighbor_offset_begin,
  BandNeighborOffsetIt const band_neighbor_offset_end)
{
  using namespace std;

  enum class DilationCell : uint8_t {
    kBackground = uint32_t{0},
    kForeground,
    kDilated
  };

  assert(LinearSize(grid_size) > size_t{0});

  if (grid_indices.empty()) {
    return vector<vector<array<int32_t, N>>>();;
  }

  // Dilation grid is padded one cell in each dimension.
  // Then transform the provided indices to the dilation grid.
  auto dilation_grid_size = grid_size;
  for_each(begin(dilation_grid_size), end(dilation_grid_size),
           [](auto& d) { d += 2; });
  auto dilation_buffer = vector<DilationCell>(
    LinearSize(dilation_grid_size), DilationCell::kBackground);
  auto dilation_grid = Grid<DilationCell, N>(
    dilation_grid_size, dilation_buffer);
  auto dilation_grid_indices = vector<array<int32_t, N>>(grid_indices.size());
  transform(
    begin(grid_indices),
    end(grid_indices),
    begin(dilation_grid_indices),
    [=](auto const& grid_index) {
      assert(Inside(grid_index, grid_size));
      auto const dilation_grid_index =
        DilationGridIndexFromDistanceGridIndex(grid_index);
      assert(Inside(dilation_grid_index, dilation_grid_size));
      return dilation_grid_index;
    });

  // Set foreground from the provided (transformed) indices.
  for (auto const& dilation_grid_index : dilation_grid_indices) {
    assert(Inside(dilation_grid_index, dilation_grid.size()) && "Precondition");
    dilation_grid.Cell(dilation_grid_index) = DilationCell::kForeground;
  }

  // Tag background cells connected to foreground as dilated.
  // We only overwrite background cells here.
  auto dilation_indices = vector<array<int32_t, N>>();
  for (auto const& dilation_grid_index : dilation_grid_indices) {
    assert(dilation_grid.Cell(dilation_grid_index) == DilationCell::kForeground);
    for (auto dilation_neighbor_offset_iter = dilation_neighbor_offset_begin;
         dilation_neighbor_offset_iter != dilation_neighbor_offset_end;
         ++dilation_neighbor_offset_iter) {
      auto neighbor_index = dilation_grid_index;
      for (auto i = size_t{0}; i < N; ++i) {
        neighbor_index[i] += (*dilation_neighbor_offset_iter)[i];
      }

      assert(Inside(neighbor_index, dilation_grid.size()) && "Precondition");
      if (dilation_grid.Cell(neighbor_index) == DilationCell::kBackground) {
        dilation_grid.Cell(neighbor_index) = DilationCell::kDilated;
        dilation_indices.push_back(neighbor_index);
      }
    }
  }
  assert(!dilation_indices.empty());

  // Get connected components of dilated cells.
  auto const dilation_bands = ConnectedComponents(
    dilation_indices,
    dilation_grid_size,
    band_neighbor_offset_begin,
    band_neighbor_offset_end);
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
template<std::size_t N> inline
std::vector<std::array<std::int32_t, N>>
NarrowBandDilationBandCells(
  std::vector<std::array<std::int32_t, N>> const& dilation_band_indices,
  Grid<uint8_t, N> const& boundary_mask_grid)
{
  using namespace std;

  if (dilation_band_indices.empty()) {
    return vector<array<int32_t, N>>();
  }

  auto narrow_band_indices = vector<array<int32_t, N>>();
  narrow_band_indices.reserve(dilation_band_indices.size());
  for (auto const& dilation_grid_index : dilation_band_indices) {
    // Since dilation bands are constructed using a vertex neighborhood,
    // not all dilation cells are face-connected to a boundary cell.
    // We add only those dilation cells that are face-connected to a
    // boundary cell, since this will be required when estimating distance
    // (i.e. solving the eikonal equation).
    auto const distance_grid_index =
      DistanceGridIndexFromDilationGridIndex(dilation_grid_index);
    assert(boundary_mask_grid.Cell(distance_grid_index) != uint8_t{1});

    // If the distance grid index is not inside the boundary mask
    // (i.e. distance) grid it cannot belong to a narrow band.
    if (Inside(distance_grid_index, boundary_mask_grid.size())) {
      // Check for boundary face-neighbors in each dimension.
      // If we find one boundary face-neighbor we are done.
      for (auto i = size_t{0}; i < N; ++i) {
        // +1
        auto neighbor_index = distance_grid_index;
        neighbor_index[i] += int32_t{1};
        if (Inside(neighbor_index, boundary_mask_grid.size()) &&
            boundary_mask_grid.Cell(neighbor_index) == uint8_t{1}) {
          narrow_band_indices.push_back(distance_grid_index);
          break;
        }
        // -1
        neighbor_index = distance_grid_index;
        neighbor_index[i] -= int32_t{1};
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


//!
//!
template<std::size_t N> inline
std::pair<std::vector<std::array<std::int32_t, N>>,
          std::vector<std::array<std::int32_t, N>>>
SignedNarrowBandIndices(
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  std::array<std::size_t, N> const& distance_grid_size)
{
  using namespace std;

  auto inside_narrow_band_indices = vector<array<int32_t, N>>();
  auto outside_narrow_band_indices = vector<array<int32_t, N>>();
  if (boundary_indices.empty()) {
    return {outside_narrow_band_indices, inside_narrow_band_indices};
  }

  // Compute connected components of boundary cells.
  auto const vtx_neighbor_offsets = VertexNeighborOffsets<N>();
  auto const connected_components = ConnectedComponents(
    boundary_indices,
    distance_grid_size,
    begin(vtx_neighbor_offsets),
    end(vtx_neighbor_offsets));
  assert(!connected_components.empty());

  // Create a mask where:
  // - boundary cells = 1
  // - non-boundary cells = 0
  auto boundary_mask_buffer =
    vector<uint8_t>(LinearSize(distance_grid_size), uint8_t{0});
  auto boundary_mask_grid =
    Grid<uint8_t, N>(distance_grid_size, boundary_mask_buffer);
  for (auto const& boundary_index : boundary_indices) {
    assert(Inside(boundary_index, boundary_mask_grid.size()) && "Precondition");
    boundary_mask_grid.Cell(boundary_index) = uint8_t{1};
  }

  // Check dilation bands of connected boundary components.
  // Dilation bands must be computed per connected component since each
  // component has a separate outer dilation band. If we were to compute
  // dilation bands for all boundary indices at once we would then need to
  // do work to figure out if these were outer or inner dilation bands.
  auto const face_neighbor_offsets = FaceNeighborOffsets<N>();
  for (auto const& connected_component : connected_components) {
    auto const dilation_bands = DilationBands(
      connected_component,
      distance_grid_size,
      begin(vtx_neighbor_offsets),
      end(vtx_neighbor_offsets),
      begin(face_neighbor_offsets),
      end(face_neighbor_offsets));
    assert(!dilation_bands.empty());

    if (dilation_bands.size() == 1) {
      // Only one dilation band means that the connected component has genus
      // zero, i.e. no holes. Thus, the dilation band must belong to
      // the outside.
      //
      // Note that the outer *dilation band* can never be empty, but the
      // *outer narrow band* can be!
      auto const& outer_dilation_band = dilation_bands.front();
      assert(!outer_dilation_band.empty());
      auto const outer_narrow_band_indices = NarrowBandDilationBandCells(
         outer_dilation_band,
         boundary_mask_grid);

      // Note that the outer narrow band is empty when the whole border of the
      // distance grid is frozen.
      outside_narrow_band_indices.insert(
        end(outside_narrow_band_indices),
        begin(outer_narrow_band_indices),
        end(outer_narrow_band_indices));
    }
    else {
      // We have more than one dilation band: one outer and possibly several
      // inner. The outer dilation band has the largest bounding box.
      // Note that when we have several dilation bands none of them can be
      // empty. The reasoning is that an empty dilation band requires the
      // whole distance grid to be frozen, in which case there cannot exist
      // an inner area.
      //
      // We compute the bounding boxes in dilation grid coordinates. This is
      // necessary since the entire outer dilation band may not be inside the
      // distance grid.
      auto dilation_band_areas = vector<pair<size_t, size_t>>();
      dilation_band_areas.reserve(dilation_bands.size());
      for (auto i = size_t{0}; i < dilation_bands.size(); ++i) {
        auto const& dilation_band = dilation_bands[i];
        assert(!dilation_band.empty());
        dilation_band_areas.push_back(
          {i, HyperVolume(BoundingBox(dilation_bands[i]))});
      }

      // Sort dilation bands by descending volume. The outer dilation band
      // is then the first element.
      sort(
        begin(dilation_band_areas),
        end(dilation_band_areas),
        [](auto const& lhs, auto const& rhs) {
          return lhs.second > rhs.second;
        });

      // Outer dilation bands of several connected components may overlap.
      // We are fine with adding an index multiple times to the outside
      // narrow band. The smallest distance will be used first and the rest
      // will be ignored. Worst-case we estimate distances for cells that
      // are not impactful.
      {
      auto const& outer_dilation_band_indices =
        dilation_bands[dilation_band_areas[0].first];
      assert(!outer_dilation_band_indices.empty());
      auto const outer_narrow_band_indices = NarrowBandDilationBandCells(
        outer_dilation_band_indices,
        boundary_mask_grid);
      assert(none_of(
               begin(outer_narrow_band_indices),
               end(outer_narrow_band_indices),
               [=](auto const& distance_grid_index) {
                 return !Inside(distance_grid_index, distance_grid_size);
               }));
      // Note that the outer narrow band is empty when the whole border of the
      // distance grid is frozen.
      outside_narrow_band_indices.insert(
        end(outside_narrow_band_indices),
        begin(outer_narrow_band_indices),
        end(outer_narrow_band_indices));
      }

      // Inner dilation bands cannot overlap.
      for (auto k = size_t{1}; k < dilation_band_areas.size(); ++k) {
        auto const& inner_dilation_band_indices =
          dilation_bands[dilation_band_areas[k].first];
        assert(!inner_dilation_band_indices.empty());
        auto const inner_narrow_band_indices = NarrowBandDilationBandCells(
          inner_dilation_band_indices,
          boundary_mask_grid);
        assert(!inner_narrow_band_indices.empty());
        assert(none_of(
                 begin(inner_narrow_band_indices),
                 end(inner_narrow_band_indices),
                 [=](auto const& distance_grid_index) {
                   return !Inside(distance_grid_index, distance_grid_size);
                 }));
        inside_narrow_band_indices.insert(
          end(inside_narrow_band_indices),
          begin(inner_narrow_band_indices),
          end(inner_narrow_band_indices));
      }
    }
  }

  return {outside_narrow_band_indices, inside_narrow_band_indices};
}


//! Returns a narrow band store containing estimated distances for the cells
//! in @a narrow_band_indices. Assumes that boundary condition distances
//! have already been set in @a distance_grid. Note that @a narrow_band_indices
//! may contain duplicates. Narrow band indices are assumed to be inside
//! @a distance_grid.
template<typename T, std::size_t N, typename E> inline
std::unique_ptr<NarrowBandStore<T, N>>
InitializedNarrowBand(
  std::vector<std::array<std::int32_t, N>> const& narrow_band_indices,
  Grid<T, N> const& distance_grid,
  E const& eikonal_solver)
{
  using namespace std;

  assert(!narrow_band_indices.empty());

  auto narrow_band =
    unique_ptr<NarrowBandStore<T, N>>(new NarrowBandStore<T, N>());

  for (auto const& narrow_band_index : narrow_band_indices) {
    assert(Inside(narrow_band_index, distance_grid.size()));
    assert(!Frozen(distance_grid.Cell(narrow_band_index)));

    narrow_band->Push({
      eikonal_solver.Solve(narrow_band_index, distance_grid),
      narrow_band_index});
  }

  return narrow_band;
}


//! Returns a narrow band store containing estimated distances for the
//! face-neighbors of @a boundary_indices.
//!
//! Precondition:
//! - Boundary condition distance have been set in @a distance_grid.
template<typename T, std::size_t N, typename E> inline
std::unique_ptr<NarrowBandStore<T, N>>
InitialUnsignedNarrowBand(
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  Grid<T, N> const& distance_grid,
  E const& eikonal_solver)
{
  using namespace std;

  static_assert(N > 0, "dimensionality cannot be zero");
  static_assert(N == E::kDimension, "mismatching eikonal solver dimension");

  assert(!boundary_indices.empty());

  auto narrow_band =
    unique_ptr<NarrowBandStore<T, N>>(new NarrowBandStore<T, N>());

  // Find non-frozen face neighbors of boundary cells, which then become
  // the initial narrow band.
  auto const kNeighborOffsets = array<int32_t, 2>{{-1, 1}};
  for (auto const& boundary_index : boundary_indices) {
    assert(Inside(boundary_index, distance_grid.size()));
    assert(Frozen(distance_grid.Cell(boundary_index)));

    for (auto i = size_t{0}; i < N; ++i) {
      for (auto const neighbor_offset : kNeighborOffsets) {
        auto neighbor_index = boundary_index;
        neighbor_index[i] += neighbor_offset;
        if (Inside(neighbor_index, distance_grid.size()) &&
            !Frozen(distance_grid.Cell(neighbor_index))) {
          // Found non-frozen face neighbor inside the distance grid
          // of a boundary index.
          narrow_band->Push({
            eikonal_solver.Solve(neighbor_index, distance_grid),
            neighbor_index});
        }
      }
    }
  }

  assert(!narrow_band->empty());

  return narrow_band;
}


//! Set boundary distances on @a distance_grid. Distances are multiplied by
//! @a multiplier (typically 1 or -1).
//!
//! Throws std::invalid_argument if:
//! - Not the same number of @a indices and @a distances, or
//! - @a indices (and @a distances) are empty, or
//! - Any index is outside the @a distance_grid, or
//! - Any duplicate in @a indices, or
//! - Any value in @a distances does not pass the @a distance_predicate test.
template <typename T, std::size_t N, typename D> inline
void SetBoundaryCondition(
  std::vector<std::array<std::int32_t, N>> const& indices,
  std::vector<T> const& distances,
  T const multiplier,
  D const distance_predicate,
  Grid<T, N>* const distance_grid)
{
  using namespace std;

  assert(distance_grid != nullptr);

  ThrowIfEmptyBoundaryIndices(indices);

  if (indices.size() != distances.size()) {
    throw invalid_argument("boundary indices/distances size mismatch");
  }

  auto const distance_grid_size = distance_grid->size();
  for (auto i = size_t{0}; i < indices.size(); ++i) {
    auto const index = indices[i];
    auto const distance = multiplier * distances[i];

    ThrowIfBoundaryIndexOutsideGrid(index, distance_grid_size);
    ThrowIfInvalidBoundaryDistance(distance_predicate(distance), distance);

    auto& distance_cell = distance_grid->Cell(index);
    ThrowIfDuplicateBoundaryIndex(Frozen(distance_cell), index);
    distance_cell = distance;
  }
}


//! Compute distances using the @a eikonal_solver for the face-neighbors of
//! the cell @a index. The distances are not written to the @a distance_grid,
//! but are instead stored in the @a narrow_band.
template <typename T, std::size_t N, typename E> inline
void UpdateNeighbors(
  std::array<std::int32_t, N> const& index,
  E const& eikonal_solver,
  Grid<T, N>* const distance_grid,
  NarrowBandStore<T, N>* const narrow_band)
{
  using namespace std;

  static_assert(N > 0, "dimensionality cannot be zero");
  static_assert(N == E::kDimension, "mismatching eikonal solver dimension");

  assert(distance_grid != nullptr);
  assert(narrow_band != nullptr);
  assert(Inside(index, distance_grid->size()));
  assert(Frozen(distance_grid->Cell(index)));

  // Update the narrow band. Check face-neighbors in all dimensions.
  auto const kNeighborOffsets = array<int32_t, 2>{{-1, 1}};
  for (auto i = size_t{0}; i < N; ++i) {
    for (auto const neighbor_offset : kNeighborOffsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;

      if (Inside(neighbor_index, distance_grid->size())) {
        // If the neighbor is not frozen compute a distance for it.
        // Note that we don't check if there is an entry for this index
        // in the narrow band already. If we happen to insert multiple
        // distances for the same index the smallest one will be frozen first
        // when marching and the larger distances will be ignored.
        auto& distance_cell = distance_grid->Cell(neighbor_index);
        if (!Frozen(distance_cell)) {
          narrow_band->Push({
            eikonal_solver.Solve(neighbor_index, *distance_grid),
            neighbor_index});
        }
      }
    }
  }
}


//! Compute distances using @a eikonal_solver for all non-frozen cells in
//! @a distance_grid that have a face-connected path to at least one of the
//! cells in @a narrow_band.
template <typename T, std::size_t N, typename E> inline
void MarchNarrowBand(
  E const& eikonal_solver,
  NarrowBandStore<T, N>* const narrow_band,
  Grid<T, N>* const distance_grid)
{
  using namespace std;

  assert(distance_grid != nullptr);
  assert(narrow_band != nullptr);

  // TODO: necessary?
  if (narrow_band->empty()) {
    throw invalid_argument("cannot march empty narrow band");
  }

  while (!narrow_band->empty()) {
    // Take smallest distance from the narrow band and freeze it.
    auto const narrow_band_cell = narrow_band->Pop();
    auto const distance = narrow_band_cell.first;
    auto const index = narrow_band_cell.second;
    assert(Frozen(distance));
    assert(Inside(index, distance_grid->size()));

    auto& distance_cell = distance_grid->Cell(index);

    // Since we allow multiple distances for the same cell index in the
    // narrow band it could be that the distance for this grid cell has
    // already been frozen. In that case just ignore subsequent values from
    // the narrow band for that grid cell and move on.
    if (!Frozen(distance_cell)) {
      distance_cell = distance;

      // Update distances for non-frozen face-neighbors of the newly
      // frozen cell.
      UpdateNeighbors(index, eikonal_solver, distance_grid, narrow_band);
    }
  }
}


//! Polynomial coefficients are equivalent to array index,
//! i.e. Sum(q[i] * x^i) = 0, for i in [0, 2], or simpler
//! q[0] + q[1] * x + q[2] * x^2 = 0.
//!
//! Returns the largest real root.
//!
//! Throws a std::runtime_error if no real roots exist.
template<typename T>
T SolveEikonalQuadratic(std::array<T, 3> const& q)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "quadratic coefficients must be floating point");

  assert(fabs(q[2]) > T(1e-9));

  auto const discriminant = q[1] * q[1] - T{4} * q[2] * q[0];
  if (discriminant < T{0}) {
    throw runtime_error("negative discriminant");
  }

  auto const root = (-q[1] + sqrt(discriminant)) / (T{2} * q[2]);
  assert(!isnan(root));

  if (root < T{0}) {
    throw runtime_error("negative distance");
  }

  return root;
}


template<typename T, std::size_t N>
T SolveEikonal(
  std::array<std::int32_t, N> const& index,
  Grid<T, N> const& distance_grid,
  T const speed,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(Inside(index, distance_grid.size()));

  auto const neighbor_offsets = array<int32_t, 2>{{-1, 1}};

  // Initialize quadratic coefficients.
  auto q = array<T, 3>{{T{-1} / Squared(speed), T{0}, T{0}}};

  // Find the smallest frozen neighbor (if any) in each dimension.
  for (auto i = size_t{0}; i < N; ++i) {
    auto neighbor_min_distance = numeric_limits<T>::max();

    // Check neighbors in both directions for this dimenion.
    for (auto const neighbor_offset : neighbor_offsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;
      if (Inside(neighbor_index, distance_grid.size())) {
        // Note that if the neighbor is not frozen it will have the default
        // distance numeric_limits<T>::max().
        auto const neighbor_distance = distance_grid.Cell(neighbor_index);
        if (neighbor_distance < neighbor_min_distance) {
          assert(Frozen(neighbor_distance));
          neighbor_min_distance = neighbor_distance;
        }
      }
    }

    // Update quadratic coefficients for the current direction.
    // If no frozen neighbor was found that dimension does not contribute
    // to the coefficients.
    if (neighbor_min_distance < numeric_limits<T>::max()) {
      auto const inv_grid_spacing_squared = InverseSquared(grid_spacing[i]);
      q[0] += Squared(neighbor_min_distance) * inv_grid_spacing_squared;
      q[1] += T{-2} * neighbor_min_distance * inv_grid_spacing_squared;
      q[2] += inv_grid_spacing_squared;
    }
  }

  return SolveEikonalQuadratic(q);
}


template<typename T, std::size_t N>
T HighAccuracySolveEikonal(
  std::array<std::int32_t, N> const& index,
  Grid<T, N> const& distance_grid,
  T const speed,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(Inside(index, distance_grid.size()));

  static auto const neighbor_offsets = array<int32_t, 2>{{-1, 1}};

  // Initialize quadratic coefficients.
  auto q = array<T, 3>{{T{-1} / Squared(speed), T{0}, T{0}}};

  // Find the smallest frozen neighbor(s) (if any) in each dimension.
  for (auto i = size_t{0}; i < N; ++i) {
    auto neighbor_min_distance = numeric_limits<T>::max();
    auto neighbor_min_distance2 = numeric_limits<T>::max();

    // Check neighbors in both directions for this dimenion.
    for (auto const neighbor_offset : neighbor_offsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;
      if (Inside(neighbor_index, distance_grid.size())) {
        auto const neighbor_distance = distance_grid.Cell(neighbor_index);
        if (neighbor_distance < neighbor_min_distance) {
          // Neighbor one step away is frozen.
          neighbor_min_distance = neighbor_distance;

          // Check if neighbor two steps away is frozen and has smaller
          // (or equal) distance than neighbor one step away.
          auto neighbor_index2 = neighbor_index;
          neighbor_index2[i] += neighbor_offset;
          if (Inside(neighbor_index2, distance_grid.size())) {
            auto const neighbor_distance2 = distance_grid.Cell(neighbor_index2);
            if (neighbor_distance2 <= neighbor_distance) {
              // Neighbor index two steps away is frozen.
              neighbor_min_distance2 = neighbor_distance2;
            }
          }
        }
      }
    }

    // Update quadratic coefficients for the current direction.
    if (neighbor_min_distance < numeric_limits<T>::max()) {
      if (neighbor_min_distance2 < numeric_limits<T>::max()) {
        // Second order coefficients.
        auto const alpha = T{9} / (T{4} * Squared(grid_spacing[i]));
        auto const t = (T{1} / T{3}) * (T{4} * neighbor_min_distance - neighbor_min_distance2);
        q[0] += Squared(t) * alpha;
        q[1] += T{-2} * t * alpha;
        q[2] += alpha;
      }
      else {
        // First order coefficients.
        auto const inv_grid_spacing_squared = InverseSquared(grid_spacing[i]);
        q[0] += Squared(neighbor_min_distance) * inv_grid_spacing_squared;
        q[1] += T{-2} * neighbor_min_distance * inv_grid_spacing_squared;
        q[2] += inv_grid_spacing_squared;
      }
    }
  }

  return SolveEikonalQuadratic(q);
}


//! Base class for Eikonal solvers.
template <typename T, std::size_t N>
class EikonalSolverBase
{
public:
  typedef T ScalarType;
  static std::size_t const kDimension = N;

protected:
  explicit EikonalSolverBase(std::array<T, N> const& grid_spacing)
    : grid_spacing_(grid_spacing)
  {
    ThrowIfInvalidGridSpacing(grid_spacing_);
  }

  std::array<T, N> grid_spacing() const
  {
    return grid_spacing_;
  }

private:
  std::array<T, N> const grid_spacing_;
};


//! Base class for Eikonal solvers with uniform speed.
template <typename T, std::size_t N>
class UniformSpeedEikonalSolverBase : public EikonalSolverBase<T, N>
{
protected:
  UniformSpeedEikonalSolverBase(
    std::array<T, N> const& grid_spacing,
    T const speed)
    : EikonalSolverBase<T, N>(grid_spacing)
    , speed_(speed)
  {
    ThrowIfZeroOrNegativeOrNanSpeed(speed_);
  }

  //! Returns the uniform speed, guaranteed to be:
  //! - Non-zero
  //! - Non-negative
  //! - Non-NaN
  T speed() const
  {
    return speed_;
  }

private:
  T const speed_;
};


//! Base class for Eikonal solvers with varying speed.
template <typename T, std::size_t N>
class VaryingSpeedEikonalSolverBase : public EikonalSolverBase<T, N>
{
protected:
  VaryingSpeedEikonalSolverBase(
    std::array<T, N> const& grid_spacing,
    std::array<std::size_t, N> const& speed_grid_size,
    std::vector<T> const& speed_buffer)
    : EikonalSolverBase(grid_spacing)
    , speed_grid_(speed_grid_size, speed_buffer)
  {
    for (auto const speed : speed_buffer) {
      ThrowIfZeroOrNegativeOrNanSpeed(speed);
    }
  }

  //! Returns the speed at @a index, guaranteed to be:
  //! - Non-zero
  //! - Non-negative
  //! - Non-NaN
  //!
  //! Throws a std::invalid_argument if @a index is outside the
  //! underlying grid.
  T Speed(std::array<std::int32_t, N> const& index) const
  {
    using namespace std;

    if (!Inside(index, speed_grid_.size())) {
      throw invalid_argument("index outside speed grid");
    }

    return speed_grid_.Cell(index);
  }

private:
  ConstGrid<T, N> const speed_grid_;
};


//!
//!
//!
template<typename T, std::size_t N, typename EikonalSolverType> inline
std::vector<T> SignedDistance(
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  std::vector<T> const& boundary_distances,
  std::vector<std::array<std::int32_t, N>> const& inside_narrow_band_indices,
  std::vector<std::array<std::int32_t, N>> const& outside_narrow_band_indices,
  EikonalSolverType const& eikonal_solver)
{
  using namespace std;

  typedef T DistanceType;

  static_assert(N > 0, "dimensionality cannot be zero");
  static_assert(N == EikonalSolverType::kDimension,
                "mismatching eikonal solver dimension");

  assert(!(inside_narrow_band_indices.empty() &&
           outside_narrow_band_indices.empty()));

  auto distance_buffer = vector<DistanceType>(
    LinearSize(grid_size), numeric_limits<DistanceType>::max());
  auto distance_grid = Grid<DistanceType, N>(grid_size, distance_buffer);
  auto distance_predicate = [](auto const d) {
    return !isnan(d) && Frozen(d);
  };

  assert(none_of(begin(distance_buffer), end(distance_buffer),
                 [](DistanceType const d) { return Frozen(d); }));

  if (!inside_narrow_band_indices.empty()) {
    // Set boundaries for marching inside (negative distance)
    SetBoundaryCondition(
      boundary_indices,
      boundary_distances,
      DistanceType{-1}, // Multiplier
      distance_predicate,
      &distance_grid);

    // Initialize inside narrow band with negated boundary distances.
    auto inside_narrow_band = InitializedNarrowBand(
      inside_narrow_band_indices,
      distance_grid,
      eikonal_solver);
    MarchNarrowBand(eikonal_solver, inside_narrow_band.get(), &distance_grid);

    // Negate all the inside distance values and set the boundary cells to have
    // their actual values. Essentially, negate everything computed so far.
    for_each(
      begin(distance_buffer),
      end(distance_buffer),
      [](auto& d) {
        d = Frozen(d) ? d * DistanceType{-1} : d;
      });
  }
  else {
    // Set boundaries for marching outside (positive distance)
    SetBoundaryCondition(
      boundary_indices,
      boundary_distances,
      DistanceType{1}, // Multiplier
      distance_predicate,
      &distance_grid);
  }

  if (!outside_narrow_band_indices.empty()) {
    // Initialize outside narrow band with "non-negated" distances.
    auto outside_narrow_band = InitializedNarrowBand(
      outside_narrow_band_indices,
      distance_grid,
      eikonal_solver);
    MarchNarrowBand(eikonal_solver, outside_narrow_band.get(), &distance_grid);
  }

  assert(all_of(begin(distance_buffer), end(distance_buffer),
                [](DistanceType const d) { return Frozen(d); }));

  return distance_buffer;

}

} // namespace detail


//! Holds parameters related to the grid and provides methods for solving
//! the eikonal equation for a single grid cell at a time using
//! the current distance grid.
template<typename T, std::size_t N>
class UniformSpeedEikonalSolver :
  public detail::UniformSpeedEikonalSolverBase<T, N>
{
public:
  explicit UniformSpeedEikonalSolver(
    std::array<T, N> const& grid_spacing,
    T const speed = T{1})
    : detail::UniformSpeedEikonalSolverBase<T, N>(grid_spacing, speed)
  {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) and states (@a state_grid) of other cells.
  T Solve(
    std::array<std::int32_t, N> const& index,
    detail::Grid<T, N> const& distance_grid) const
  {
    return detail::SolveEikonal(
      index,
      distance_grid,
      speed(),
      grid_spacing());
  }
};


template <typename T, std::size_t N>
class HighAccuracyUniformSpeedEikonalSolver :
  public detail::UniformSpeedEikonalSolverBase<T, N>
{
public:
  explicit HighAccuracyUniformSpeedEikonalSolver(
    std::array<T, N> const& grid_spacing,
    T const speed = T{1})
    : detail::UniformSpeedEikonalSolverBase<T, N>(grid_spacing, speed)
  {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) and states (@a state_grid) of other cells.
  T Solve(
    std::array<std::int32_t, N> const& index,
    detail::Grid<T, N> const& distance_grid) const
  {
    return detail::HighAccuracySolveEikonal(
      index,
      distance_grid,
      speed(),
      grid_spacing());
  }
};


template <typename T, std::size_t N>
class VaryingSpeedEikonalSolver :
  public detail::VaryingSpeedEikonalSolverBase<T, N>
{
public:
  VaryingSpeedEikonalSolver(
    std::array<T, N> const& grid_spacing,
    std::array<std::size_t, N> const& speed_grid_size,
    std::vector<T> const& speed_buffer)
    : detail::VaryingSpeedEikonalSolverBase<T, N>(
        grid_spacing, speed_grid_size, speed_buffer)
  {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) and states (@a state_grid) of other cells.
  T Solve(
    std::array<std::int32_t, N> const& index,
    detail::Grid<T, N> const& distance_grid) const
  {
    return detail::SolveEikonal(
      index,
      distance_grid,
      Speed(index),
      grid_spacing());
  }
};


template <typename T, std::size_t N>
class HighAccuracyVaryingSpeedEikonalSolver :
  public detail::VaryingSpeedEikonalSolverBase<T, N>
{
public:
  HighAccuracyVaryingSpeedEikonalSolver(
    std::array<T, N> const& grid_spacing,
    std::array<std::size_t, N> const& speed_grid_size,
    std::vector<T> const& speed_buffer)
    : detail::VaryingSpeedEikonalSolverBase<T, N>(
        grid_spacing, speed_grid_size, speed_buffer)
  {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) and states (@a state_grid) of other cells.
  T Solve(
    std::array<std::int32_t, N> const& index,
    detail::Grid<T, N> const& distance_grid) const
  {
    return detail::HighAccuracySolveEikonal(
      index,
      distance_grid,
      Speed(index),
      grid_spacing());
  }
};


//! TODO - example usage!
template<typename T, std::size_t N, typename EikonalSolverType> inline
std::vector<T> UnsignedDistance(
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  std::vector<T> const& boundary_distances,
  EikonalSolverType const& eikonal_solver)
{
  using namespace std;
  using namespace detail;

  typedef T DistanceType;

  auto distance_buffer = vector<DistanceType>(
    LinearSize(grid_size), numeric_limits<DistanceType>::max());
  auto distance_grid = Grid<DistanceType, N>(grid_size, distance_buffer);
  auto distance_predicate = [](auto const d) {
    return !isnan(d) && Frozen(d) && d >= T{0};
  };

  assert(none_of(begin(distance_buffer), end(distance_buffer),
                 [](DistanceType const d) { return Frozen(d); }));

  SetBoundaryCondition(
    boundary_indices,
    boundary_distances,
    DistanceType{1}, // Distance multiplier.
    distance_predicate,
    &distance_grid);

  auto narrow_band = InitialUnsignedNarrowBand(
    boundary_indices, distance_grid, eikonal_solver);

  MarchNarrowBand(eikonal_solver, narrow_band.get(), &distance_grid);

  assert(all_of(begin(distance_buffer), end(distance_buffer),
                [](DistanceType const d) { return Frozen(d); }));

  return distance_buffer;
}


//! Compute the signed distance on a grid.
//!
//! Input:
//!   grid_size          - Number of grid cells in each dimension.
//!   boundary_indices   - Integer coordinates of cells with provided distances.
//!   boundary_distances - Signed distances assigned to boundary cells.
//!
//! Preconditions:
//!   - grid_size may not have a zero element.
//!   - dx must have all positive elements.
//!   - speed?
//!   - frozen_indices, frozen_distances and normals must have the same size.
//!   - frozen_indices must all be within size.
//!
//! TODO - example usage!
template<typename T, std::size_t N, typename EikonalSolverType> inline
std::vector<T> SignedDistance(
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  std::vector<T> const& boundary_distances,
  EikonalSolverType const& eikonal_solver)
{
  detail::ThrowIfZeroElementInSize(grid_size);
  detail::ThrowIfEmptyBoundaryIndices(boundary_indices);
  detail::ThrowIfBoundaryIndexOutsideGrid(boundary_indices, grid_size);

  auto narrow_band_indices = detail::SignedNarrowBandIndices(
    boundary_indices,
    grid_size);
  auto outside_narrow_band_indices = narrow_band_indices.first;
  auto inside_narrow_band_indices = narrow_band_indices.second;

  return detail::SignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    inside_narrow_band_indices,
    outside_narrow_band_indices,
    eikonal_solver);
}

} // namespace fast_marching_method
} // namespace thinks

#endif // THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED
