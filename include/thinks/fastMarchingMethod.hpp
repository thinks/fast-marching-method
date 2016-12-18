#ifndef THINKS_FMM_FASTMARCHINGMETHOD_HPP_INCLUDED
#define THINKS_FMM_FASTMARCHINGMETHOD_HPP_INCLUDED

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

//! Returns the product of the elements in array @a a.
//! Note: Not checking for integer overflow here!
template<std::size_t N> inline
std::size_t LinearSize(std::array<std::size_t, N> const& a)
{
  using namespace std;

  return accumulate(begin(a), end(a), size_t{1}, multiplies<size_t>());
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


//! Throws a std::invalid_argument if one or more of the elements in
//! @a grid_size is zero.
template<std::size_t N> inline
void ThrowIfZeroElementInGridSize(std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  if (find_if(begin(grid_size), end(grid_size),
              [](auto const x) { return x == size_t{0}; }) != end(grid_size)) {
    auto ss = stringstream();
    ss << "invalid grid size: " << ToString(grid_size);
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
    ss << "boundary index outside grid: " << ToString(boundary_index);
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
    ThrowIfZeroElementInGridSize(size);
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
    ThrowIfZeroElementInGridSize(size);
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


template<std::size_t N, typename NeighborOffsetIt> inline
std::vector<std::vector<std::array<std::int32_t, N>>>
ConnectedComponents(
  std::vector<std::array<std::int32_t, N>> const& indices,
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

  ThrowIfZeroElementInGridSize(grid_size);

  auto connected_components = vector<vector<array<int32_t, N>>>();

  if (indices.empty()) {
    return connected_components;
  }

  auto label_buffer =
    vector<LabelCell>(LinearSize(grid_size), LabelCell::kBackground);
  auto label_grid = Grid<LabelCell, N>(grid_size, label_buffer);

  for (auto const& index : indices) {
    ThrowIfBoundaryIndexOutsideGrid(index, label_grid.size());
    // Note: We don't check for duplicate indices here.
    label_grid.Cell(index) = LabelCell::kForeground;
  }

  for (auto const& index : indices) {
    assert(Inside(index, label_grid.size()));
    assert(label_grid.Cell(index) == LabelCell::kForeground ||
           label_grid.Cell(index) == LabelCell::kLabelled);

    if (label_grid.Cell(index) == LabelCell::kForeground) {
      // This index has not already been labelled.
      // Start a new component.
      auto component = vector<array<int32_t, N>>();
      auto neighbor_indices = stack<array<int32_t, N>>();
      label_grid.Cell(index) = LabelCell::kLabelled;
      component.push_back(index);
      neighbor_indices.push(index);

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
            // Mark neighbor as labelled, store in component and add
            // to list of indices whose neighbors we should check.
            label_grid.Cell(neighbor_index) = LabelCell::kLabelled;
            component.push_back(neighbor_index);
            neighbor_indices.push(neighbor_index);
          }
        }
      }
      connected_components.push_back(component);
    }
  }

  return connected_components;
}


//!
template<std::size_t N> inline
std::vector<std::vector<std::array<std::int32_t, N>>>
DilationBands(
  std::vector<std::array<std::int32_t, N>> const& grid_indices,
  std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  enum class DilationCell : uint8_t {
    kBackground = uint32_t{0},
    kForeground,
    kDilated
  };

  assert(LinearSize(grid_size) > size_t{0});

  auto dilation_bands = vector<vector<array<int32_t, N>>>();
  if (grid_indices.empty()) {
    return dilation_bands;
  }

  // Dilation grid is padded one cell in each dimension (positive and negative).
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
      auto dilation_grid_index = grid_index;
      for_each(begin(dilation_grid_index), end(dilation_grid_index),
               [](auto& d) { d += int32_t{1}; });
      assert(Inside(dilation_grid_index, dilation_grid_size));
      return dilation_grid_index;
    });

  // Set foreground from provided indices.
  for (auto const& dilation_grid_index : dilation_grid_indices) {
    assert(Inside(dilation_grid_index, dilation_grid.size()));
    dilation_grid.Cell(dilation_grid_index) = DilationCell::kForeground;
  }

  // Tag background cells connected to foreground as dilated.
  // We only overwrite background cells here.
  auto const dilation_neighbor_offsets = VertexNeighborOffsets<N>();
  auto dilation_indices = vector<array<int32_t, N>>();
  for (auto const& dilation_grid_index : dilation_grid_indices) {
    assert(dilation_grid.Cell(dilation_grid_index) == DilationCell::kForeground);
    for (auto const dilation_neighbor_offset : dilation_neighbor_offsets) {
      auto neighbor_index = dilation_grid_index;
      for (auto i = size_t{0}; i < N; ++i) {
        neighbor_index[i] += dilation_neighbor_offset[i];
      }

      if (dilation_grid.Cell(neighbor_index) == DilationCell::kBackground) {
        dilation_grid.Cell(neighbor_index) = DilationCell::kDilated;
        dilation_indices.push_back(neighbor_index);
      }
    }
  }
  assert(!dilation_indices.empty());

  // Get connected components for dilated cells.
  auto const face_neighbor_offsets = FaceNeighborOffsets<N>();
  auto const connected_dilation_components = ConnectedComponents(
    dilation_indices,
    dilation_grid_size,
    begin(face_neighbor_offsets),
    end(face_neighbor_offsets));

  // A dilation band is defined as a (face) connected components of
  // dilated cells.
  for (auto const& dilation_component : connected_dilation_components) {
    auto dilation_band = vector<array<int32_t, N>>();
    for (auto const& dilation_index : dilation_component) {
      // Transform indices from dilation grid to original grid.
      // The transformed index may be outside the original grid.
      auto grid_index = dilation_index;
      for_each(begin(grid_index), end(grid_index),
               [](auto& d) { d -= int32_t{1}; });
      if (Inside(grid_index, grid_size)) {
        // Since dilation bands are constructed using a vertex neighborhood,
        // not all dilation cells are face-connected to a foreground cell.
        // We add only those dilation cells that are face-connected to a
        // foreground cell.
        auto found_foreground_face_neighbor = false;
        for (auto i = size_t{0}; i < N; ++i) {
          auto neighbor_index = dilation_index;
          neighbor_index[i] += int32_t{1};
          if (dilation_grid.Cell(neighbor_index) == DilationCell::kForeground) {
            found_foreground_face_neighbor = true;
            break;
          }
          neighbor_index = dilation_index;
          neighbor_index[i] -= int32_t{1};
          if (dilation_grid.Cell(neighbor_index) == DilationCell::kForeground) {
            found_foreground_face_neighbor = true;
            break;
          }
        }

        if (found_foreground_face_neighbor) {
          dilation_band.push_back(grid_index);
        }
      }
    }

    // NOTE!
    // We add the dilation band here even though it may be empty. This can
    // happen if all cells for a dilated component are outside the original
    // grid.
    dilation_bands.push_back(dilation_band);
  }

  return dilation_bands;
}


//! Returns an array of pairs, where each element is the min/max index
//! coordinates in the corresponding dimension.
//!
//! Throws std::invalid_argument exception if @a indices is empty.
template<std::size_t N> inline
std::array<std::pair<std::int32_t, std::int32_t>, N> BoundingBox(
  std::vector<std::array<std::int32_t, N>> const& indices)
{
  using namespace std;

  if (indices.empty()) {
    throw invalid_argument("cannot compute bounding box from empty indices");
  }

  auto bbox = array<pair<int32_t, int32_t>, N>();

  // Initialize bounding box in all dimensions.
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


//! Returns true if the (distance) value @a d indicates that a distance cell is
//! frozen, otherwise false.
template <typename T> inline
bool frozen(T const d)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(!isnan(d));
  return -numeric_limits<T>::max() < d && d < numeric_limits<T>::max();
}


//!
template<std::size_t N> inline
std::pair<std::vector<std::array<int32_t, N>>,
          std::vector<std::array<int32_t, N>>>
SignedNarrowBandIndices(
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  std::array<std::size_t, N> const& distance_grid_size)
{
  using namespace std;

  ThrowIfEmptyBoundaryIndices(boundary_indices);

  // Check dilation bands of connected boundary components.
  auto const cc_neighbor_offsets = VertexNeighborOffsets<N>();
  auto const connected_components = ConnectedComponents(
    boundary_indices,
    distance_grid_size,
    begin(cc_neighbor_offsets),
    end(cc_neighbor_offsets));
  assert(!connected_components.empty());

  auto inside_narrow_band_indices = vector<array<int32_t, N>>();
  auto outside_narrow_band_indices = vector<array<int32_t, N>>();

  for (auto const& connected_component : connected_components) {
    auto const dilation_bands = DilationBands(
      connected_component, distance_grid_size);

    assert(!dilation_bands.empty());
    if (dilation_bands.size() == 1) {
      //throw invalid_argument("open boundary");
      auto const& outer_dilation_band = dilation_bands.front();
      for (auto const& outer_dilation_index : outer_dilation_band) {
        assert(Inside(outer_dilation_index, distance_grid_size));
        outside_narrow_band_indices.push_back(outer_dilation_index);
      }
    }
    else {
      // We have more than one dilation band: one outer and possibly several
      // inner. The outer dilation band has the largest bounding box. We
      // ignore empty dilation bands here, they are used only to get the
      // corrent number of dilation bands.
      auto dilation_band_areas = vector<pair<size_t, size_t>>();
      for (auto i = size_t{0}; i < dilation_bands.size(); ++i) {
        auto const& dilation_band = dilation_bands[i];
        if (!dilation_band.empty()) {
          dilation_band_areas.push_back(
            {i, HyperVolume(BoundingBox(dilation_bands[i]))});
        }
      }
      assert(!dilation_band_areas.empty());

      // Sort dilation bands by descending volume. The outer dilation band
      // is then the first element.
      sort(
        begin(dilation_band_areas),
        end(dilation_band_areas),
        [](auto const& lhs, auto const& rhs) {
          return lhs.second > rhs.second;
        });
      auto const& outer_dilation_band =
        dilation_bands[dilation_band_areas[0].first];
      assert(!outer_dilation_band.empty());

      // Outer dilation bands of several connected components may overlap.
      // We are fine with adding an index multiple times to the outside
      // narrow band. The smallest distance will be used first and the rest
      // will be ignored.
      for (auto const& outer_dilation_index : outer_dilation_band) {
        assert(Inside(outer_dilation_index, distance_grid_size));
        outside_narrow_band_indices.push_back(outer_dilation_index);
      }

      // Inner dilation bands cannot overlap.
      for (auto k = size_t{1}; k < dilation_bands.size(); ++k) {
        auto const& inner_dilation_band =
          dilation_bands[dilation_band_areas[k].first];
        for (auto const& inner_dilation_index : inner_dilation_band) {
          assert(Inside(inner_dilation_index, distance_grid_size));
          inside_narrow_band_indices.push_back(inner_dilation_index);
        }
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
    assert(!frozen(distance_grid.Cell(narrow_band_index)));

    narrow_band->Push({
      eikonal_solver.Solve(narrow_band_index, distance_grid),
      narrow_band_index});
  }

  return narrow_band;
}


//! Returns a narrow band store containing estimated distances for the face
//! neighbors of @a boundary_indices. Assumes that boundary condition distances
//! have already been set in @a distance_grid.
template<typename T, std::size_t N, typename E> inline
std::unique_ptr<NarrowBandStore<T, N>>
InitialUnsignedNarrowBand(
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  Grid<T, N> const& distance_grid,
  E const& eikonal_solver)
{
  using namespace std;

  assert(!boundary_indices.empty());

  auto narrow_band =
    unique_ptr<NarrowBandStore<T, N>>(new NarrowBandStore<T, N>());

  // Find non-frozen face neighbors of boundary cells, which then become
  // the initial narrow band.
  auto const kNeighborOffsets = array<int32_t, 2>{{-1, 1}};
  for (auto const& boundary_index : boundary_indices) {
    assert(Inside(boundary_index, distance_grid.size()));
    assert(frozen(distance_grid.Cell(boundary_index)));

    for (auto i = size_t{0}; i < N; ++i) {
      for (auto const neighbor_offset : kNeighborOffsets) {
        auto neighbor_index = boundary_index;
        neighbor_index[i] += neighbor_offset;
        if (Inside(neighbor_index, distance_grid.size()) &&
            !frozen(distance_grid.Cell(neighbor_index))) {
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
//! - Any value in @a distances does not pass the @a distance_predicate test, or
//! - The whole grid is frozen.
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

    if (!distance_predicate(distance)) {
      auto ss = stringstream();
      ss << "invalid boundary distance: " << distance;
      throw invalid_argument(ss.str());
    }

    auto& distance_cell = distance_grid->Cell(index);
    if (frozen(distance_cell)) {
      auto ss = stringstream();
      ss << "duplicate boundary index: " << ToString(index);
      throw invalid_argument(ss.str());
    }
    distance_cell = distance;
  }

  // Here we know that all boundary indices are unique and inside the grid.
  if (indices.size() == LinearSize(distance_grid->size())) {
    throw invalid_argument("whole grid is boundary");
  }
}


//! Compute distances using the @a eikonal_solver for the face neighbors of
//! @a index. These are not written to the @a distance_grid, but are instead
//! stored in the @a narrow_band.
template <typename T, std::size_t N, typename E> inline
void UpdateNeighbors(
  std::array<std::int32_t, N> const& index,
  E const& eikonal_solver,
  Grid<T, N>* const distance_grid,
  NarrowBandStore<T, N>* const narrow_band)
{
  using namespace std;

  assert(distance_grid != nullptr);
  assert(narrow_band != nullptr);
  assert(Inside(index, distance_grid->size()));

  // Update the narrow band.
  const auto kNeighborOffsets = array<int32_t, 2>{{-1, 1}};
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
        if (!frozen(distance_cell)) {
          auto const distance = eikonal_solver.Solve(
            neighbor_index,
            *distance_grid);
          narrow_band->Push({distance, neighbor_index});
        }
      }
    }
  }
}


//! Compute distances for all cells in @a distance_grid. Starting from the
//! initial indices in @a narrow_band, freeze the smallest distance and update
//! the (non-frozen) neighbor distances for that cell and add the neighbors
//! to the @a narrow_band. Repeat this process until there are no non-frozen
//! cells in @a distance_grid.
template <typename T, std::size_t N, typename E> inline
void MarchNarrowBand(
  E const& eikonal_solver,
  NarrowBandStore<T, N>* const narrow_band,
  Grid<T, N>* const distance_grid)
{
  using namespace std;

  assert(distance_grid != nullptr);
  assert(narrow_band != nullptr);

  while (!narrow_band->empty()) {
    // Take smallest distance from the narrow band and freeze it.
    auto const narrow_band_cell = narrow_band->Pop();
    auto const distance = narrow_band_cell.first;
    auto const index = narrow_band_cell.second;
    assert(frozen(distance));
    assert(Inside(index, distance_grid->size()));

    auto& distance_cell = distance_grid->Cell(index);

    // Since we allow multiple distances for the same index in the narrow band
    // it could be that the distance for this grid cell has already been
    // frozen. In that case just ignore subsequent values from the narrow band
    // for that grid cell and move on.
    if (!frozen(distance_cell)) {
      distance_cell = distance;
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
          assert(frozen(neighbor_distance));
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

  static_assert(N > 0, "number of dimensions must be > 0");
  static_assert(N == EikonalSolverType::kDimension,
                "mismatching eikonal solver dimension");

  auto distance_buffer = vector<DistanceType>(
    LinearSize(grid_size), numeric_limits<DistanceType>::max());
  auto distance_grid = Grid<DistanceType, N>(grid_size, distance_buffer);
  auto distance_predicate = [](auto const d) {
    return !isnan(d) && frozen(d);
  };

  assert(none_of(begin(distance_buffer), end(distance_buffer),
                 [](DistanceType const d) { return frozen(d); }));

  if (!inside_narrow_band_indices.empty()) {
    // Set boundaries for marching inside (negative distance)
    SetBoundaryCondition(
      boundary_indices,
      boundary_distances,
      DistanceType{-1}, // Multiplier
      distance_predicate,
      &distance_grid);

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
        d = frozen(d) ? d * DistanceType{-1} : d;
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
    auto outside_narrow_band = InitializedNarrowBand(
      outside_narrow_band_indices,
      distance_grid,
      eikonal_solver);
    MarchNarrowBand(eikonal_solver, outside_narrow_band.get(), &distance_grid);
  }

  assert(all_of(begin(distance_buffer), end(distance_buffer),
                [](DistanceType const d) { return frozen(d); }));

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

  static_assert(N > 0, "number of dimensions must be > 0");
  static_assert(N == EikonalSolverType::kDimension,
                "mismatching eikonal solver dimension");

  auto distance_buffer = vector<DistanceType>(
    LinearSize(grid_size), numeric_limits<DistanceType>::max());
  auto distance_grid = Grid<DistanceType, N>(grid_size, distance_buffer);
  auto distance_predicate = [](auto const d) {
    return !isnan(d) && frozen(d) && d >= T{0};
  };

  assert(none_of(begin(distance_buffer), end(distance_buffer),
                 [](DistanceType const d) { return frozen(d); }));

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
                [](DistanceType const d) { return frozen(d); }));

  return distance_buffer;
}


//! Compute the signed distance on a grid.
//!
//! Input:
//!   grid_size        - Number of grid cells in each dimension.
//!   dx               - Grid cell physical size in each dimension.
//!   speed            - Interface speed, when set to one gives
//!                      Euclidean distance. Must be positive.
//!   frozen_indices   - Integer coordinates of cells with given distances.
//!   frozen_distances - Signed distances assigned to frozen cells.
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

#endif // THINKS_FMM_FASTMARCHINGMETHOD_HPP_INCLUDED
