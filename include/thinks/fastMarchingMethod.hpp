#ifndef THINKS_FMM_FASTMARCHINGMETHOD_HPP_INCLUDED
#define THINKS_FMM_FASTMARCHINGMETHOD_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <exception>
#include <future>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stack>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>


namespace thinks {
namespace fmm {
namespace detail {

//! Returns the product of the elements in array @a a.
//! Note: Not checking for integer overflow here!
template<std::size_t N> inline
std::size_t LinearSize(std::array<std::size_t, N> const& a)
{
  using namespace std;

  return accumulate(begin(a), end(a), size_t{1}, multiplies<size_t>());
}


template<typename T> inline constexpr
T Squared(T const x)
{
  return x * x;
}


template<typename T> inline constexpr
T InverseSquared(T const x)
{
  using namespace std;

  static_assert(is_floating_point<T>::value, "value must be floating point");

  return T{1} / Squared(x);
}

template<typename T, std::size_t N> inline
std::array<T, N> InverseSquared(std::array<T, N> const& a)
{
  using namespace std;

  auto r = array<T, N>();
  transform(begin(a), end(a), begin(r),
            [](T const x) { return InverseSquared(x); });
  return r;
}


#if 0
template<typename T, std::size_t N> inline
T SquaredMagnitude(std::array<T, N> const& v)
{
  using namespace std;

  auto sm = T{0};
  for (auto i = size_t{0}; i < N; ++i) {
    sm += v[i] * v[i];
  }
  return sm;
}
#endif


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


template<std::size_t N> inline
void ThrowIfInvalidGridSize(std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  if (find_if(begin(grid_size), end(grid_size),
              [](auto const x) { return x == size_t{0}; }) != end(grid_size)) {
    auto ss = stringstream();
    ss << "invalid grid size: " << ToString(grid_size);
    throw runtime_error(ss.str());
  }
}


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
    throw runtime_error(ss.str());
  }
}


template<typename T, std::size_t N> inline
void ThrowIfInvalidGridSpacing(std::array<T, N> const& grid_spacing)
{
  using namespace std;

  if (find_if(begin(grid_spacing), end(grid_spacing),
              [](auto const x) { return isnan(x) || x <= T(0); }) != end(grid_spacing)) {
    auto ss = stringstream();
    ss << "invalid grid spacing: " << ToString(grid_spacing);
    throw runtime_error(ss.str());
  }
}


template<typename T> inline
void ThrowIfInvalidSpeed(T const speed)
{
  using namespace std;

  if (isnan(speed) || speed <= T{0}) {
    auto ss = stringstream();
    ss << "invalid speed: " << speed;
    throw runtime_error(ss.str());
  }
}


template<typename T> inline
void ThrowIfInvalidSpeedBuffer(std::vector<T> const& speed_buffer)
{
  for (auto const speed : speed_buffer) {
    ThrowIfInvalidSpeed(speed);
  }
}


template<std::size_t N>
std::array<std::size_t, N - 1> GridStrides(
  std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  std::array<std::size_t, N - 1> strides;
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


//! Allows accessing a linear array as if it where a higher dimensional
//! object.
template<typename T, std::size_t N>
class Grid
{
public:
  typedef T CellType;
  typedef std::array<std::size_t, N> SizeType;
  typedef std::array<std::int32_t, N> IndexType;

  //! Construct a grid from a given @a size and @a cell_buffer.
  //!
  //! Throws a std::runtime_error if:
  //! - size evaluates to a zero linear size, i.e. if any of the elements are zero.
  //! - the linear size of @a size is not equal to the @a cell_buffer size.
  Grid(SizeType const& size, std::vector<T>& cell_buffer)
    : size_(size)
    , strides_(GridStrides(size))
    , cells_(&cell_buffer.front())
  {
    ThrowIfInvalidGridSize(size);
    ThrowIfInvalidCellBufferSize(size, cell_buffer.size());
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
  CellType* const cells_;
};


//! Allows const accessing a linear array as if it where a higher
//! dimensional object.
template<typename T, std::size_t N>
class ConstGrid
{
public:
  typedef T CellType;
  typedef std::array<std::size_t, N> SizeType;
  typedef std::array<std::int32_t, N> IndexType;

  //! Construct a grid from a given @a size and @a cell_buffer.
  //!
  //! Throws a std::runtime_error if:
  //! - size evaluates to a zero linear size, i.e. if any of the elements are zero.
  //! - the linear size of @a size is not equal to the @a cell_buffer size.
  ConstGrid(SizeType const& size, std::vector<T> const& cell_buffer)
    : size_(size)
    , strides_(GridStrides(size))
    , cells_(&cell_buffer.front())
  {
    ThrowIfInvalidGridSize(size);
    ThrowIfInvalidCellBufferSize(size, cell_buffer.size());
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
  CellType const* const cells_;
};


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
    assert(values_.empty() == index_to_pos_.empty());
    return values_.empty();
  }

  // O(log N)
  ValueType Pop()
  {
    if (empty()) {
      throw std::runtime_error("cannot pop empty narrow band store");
    }

    // Grab the top of the heap and use as return value below.
    auto const value = *values_.begin();

    // Place value from leaf level on top.
    Swap_(0, values_.size() - 1);
    index_to_pos_.erase(value.second); // ~O(1), depends on hashing.
    values_.pop_back(); // O(1)
    assert(values_.size() == index_to_pos_.size());

    // Sift the new top value downwards to restore heap constraints.
    if (!empty()) {
      SiftDown_(0);
    }

    return value;
  }

  void Insert(ValueType const& value)
  {
    using namespace std;

    if (index_to_pos_.find(value.second) != index_to_pos_.end()) {
      throw runtime_error("cannot insert existing index in narrow band store");
    }

    // Insert value at leaf level and sift it upwards.
    auto const pos = values_.size();
    values_.push_back(value);
    index_to_pos_.insert({value.second, pos});
    SiftUp_(pos);
  }

  void IncreaseDistance(
    IndexType const& index,
    DistanceType const new_distance)
  {
    using namespace std;

    auto const pos_iter = index_to_pos_.find(index);
    if (pos_iter == index_to_pos_.end()) {
      throw runtime_error("index not found");
    }

    auto& value = values_[pos_iter->second];
    if (new_distance <= value.first) {
      throw runtime_error("new distance must be greater than existing distance");
    }

    value.first = new_distance;
    SiftDown_(pos_iter->second);
  }

  void DecreaseDistance(
    IndexType const& index,
    DistanceType const new_distance)
  {
    using namespace std;

    auto const pos_iter = index_to_pos_.find(index);
    if (pos_iter == index_to_pos_.end()) {
      throw runtime_error("index not found");
    }

    auto& value = values_[pos_iter->second];
    if (new_distance >= value.first) {
      throw runtime_error("new distance must be less than existing distance");
    }

    value.first = new_distance;
    SiftUp_(pos_iter->second);
  }

private:
  typedef typename std::vector<ValueType>::size_type SizeType_;

  void SiftUp_(SizeType_ const pos)
  {
    assert(pos < values_.size());
    if (pos == 0) {
      return; // Reached the top of the heap.
    }

    // Swap "upwards" (i.e. with parent) while parent value is larger.
    auto const parent_pos = parent_pos_(pos);
    assert(parent_pos < values_.size());
    auto& pos_value = values_[pos].first;
    auto& parent_value = values_[parent_pos].first;
    if (pos_value < parent_value) {
      Swap_(pos, parent_pos);
      SiftUp_(parent_pos); // Recursive!
    }
  }

  void SiftDown_(SizeType_ const pos)
  {
    assert(pos < values_.size());
    auto const left_pos = left_child_pos_(pos);
    auto const right_pos = right_child_pos_(pos);
    auto const max_pos = values_.size() - 1;

    assert(left_pos < right_pos);
    if (left_pos > max_pos) {
      // Pos is a leaf since left child is outside array,
      // and right child is even further outside.
      return;
    }

    // Check distance values of left and right children.
    auto min_pos = pos;
    auto min_value = values_[min_pos].first;

    auto const left_value = values_[left_pos].first;
    if (left_value < min_value) {
      min_pos = left_pos;
      min_value = values_[min_pos].first;
    }

    if (right_pos <= max_pos) {
      auto const right_value = values_[right_pos].first;
      if (right_value < min_value) {
        min_pos = right_pos;
      }
    }

    // Swap with the child that has the smaller distance value,
    // if any of the child distance values is smaller than the current distance.
    if (min_pos != pos) {
      Swap_(min_pos, pos);
      SiftDown_(min_pos); // Recursive!
    }
  }

  static SizeType_ parent_pos_(SizeType_ const child_pos)
  {
    assert(child_pos > 0);
    return (child_pos - 1) / 2;
  }

  static SizeType_ left_child_pos_(SizeType_ const parent_pos)
  {
    return 2 * parent_pos + 1;
  }

  static SizeType_ right_child_pos_(SizeType_ const parent_pos)
  {
    return 2 * parent_pos + 2;
  }

  void Swap_(SizeType_ const pos0, SizeType_ const pos1)
  {
    assert(pos0 < values_.size());
    assert(pos1 < values_.size());
    auto const iter0 = index_to_pos_.find(values_[pos0].second);
    auto const iter1 = index_to_pos_.find(values_[pos1].second);
    assert(iter0 != index_to_pos_.end());
    assert(iter1 != index_to_pos_.end());
    iter0->second = pos1;
    iter1->second = pos0;
    std::swap(values_[pos0], values_[pos1]);
  }

  template<typename V, typename H> static inline
  void hash_combine_(V const& v, H const& hasher, std::size_t& seed)
  {
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  struct HashType_
  {
    typedef IndexType ArgumentType;
    typedef size_t ResultType;

    ResultType operator()(ArgumentType const& a) const
    {
      using namespace std;

      typedef typename ArgumentType::value_type ValueType;

      hash<ValueType> hasher;
      auto seed = size_t{0};
      for (auto i = size_t{0}; i < N; ++i) {
        hash_combine_(a[i], hasher, seed);
      }
      return seed;
    }
  };

  struct EqualType_
  {
    typedef bool ResultType;
    typedef IndexType FirstArgumentType;
    typedef IndexType SecondArgumentType;

    ResultType operator()(FirstArgumentType const& lhs,
                          SecondArgumentType const& rhs) const
    {
      for (auto i = size_t{0}; i < N; ++i) {
        if (lhs[i] != rhs[i]) {
          return false;
        }
      }
      return true;
    }
  };

  std::vector<ValueType> values_;
  std::unordered_map<IndexType, SizeType_, HashType_, EqualType_> index_to_pos_;
};


enum class MarchingCellState
{
  kFar = 0,
  kNarrowBand,
  kFrozen
};





//! Returns base^exponent as a compile-time constant.
constexpr std::size_t static_pow(std::size_t const base, std::size_t const exponent)
{
  using namespace std;

  // NOTE: Cannot use loops in constexpr functions in C++11, have to use
  // recursion here.
  return exponent == size_t{0} ? size_t{1} : base * static_pow(base, exponent - 1);
}


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


template<std::size_t N> inline
std::array<std::array<std::int32_t, N>, static_pow(3, N) - 1> VertexNeighborOffsets()
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
void ConnectedComponents(
  std::vector<std::array<std::int32_t, N>> const& indices,
  std::array<std::size_t, N> const& grid_size,
  NeighborOffsetIt const neighbor_offset_begin,
  NeighborOffsetIt const neighbor_offset_end,
  std::vector<std::vector<std::array<std::int32_t, N>>>* connected_components)
{
  using namespace std;

  assert(LinearSize(grid_size) > size_t{0});
  assert(connected_components != nullptr);

  connected_components->clear();
  if (indices.empty()) {
    return;
  }

  enum class LabelCell : uint8_t {
    kBackground = uint32_t{0},
    kForeground,
    kLabelled
  };

  auto label_buffer =
    vector<LabelCell>(LinearSize(grid_size), LabelCell::kBackground);
  auto label_grid = Grid<LabelCell, N>(grid_size, label_buffer.front());

  for (auto const& index : indices) {
    assert(Inside(index, label_grid.size()));
    label_grid.Cell(index) = LabelCell::kForeground;
  }

  for (auto const& index : indices) {
    assert(Inside(index, label_grid.size()));
    assert(label_grid.Cell(index) == LabelCell::kForeground ||
           label_grid.Cell(index) == LabelCell::kLabelled);
    // Check if this index has been labelled already.
    if (label_grid.Cell(index) == LabelCell::kForeground) {
      // Start a new component.
      label_grid.Cell(index) = LabelCell::kLabelled;
      auto component = vector<array<int32_t, N>>();
      component.push_back(index);
      auto neighbor_indices = stack<array<int32_t, N>>();
      neighbor_indices.push(index);

      // Flood-fill current label.
      while (!neighbor_indices.empty()) {
        auto const top_neighbor_index = neighbor_indices.top();
        neighbor_indices.pop();
        for (auto neighbor_offset_iter = neighbor_offset_begin;
             neighbor_offset_iter != neighbor_offset_end; ++neighbor_offset_iter) {
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
      connected_components->push_back(component);
    }
  }
}


template<std::size_t N> inline
void DilationBands(
  std::vector<std::array<std::int32_t, N>> const& indices,
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::vector<std::array<std::int32_t, N>>>* dilation_bands)
{
  using namespace std;

  assert(LinearSize(grid_size) > size_t{0});
  assert(dilation_bands != nullptr);

  dilation_bands->clear();
  if (indices.empty()) {
    return;
  }

  enum class DilationCell : uint8_t {
    kBackground = uint32_t{0},
    kForeground,
    kDilated
  };

  // Dilation grid is padded one cell in each dimension (positive and negative).
  auto dilation_grid_size = grid_size;
  for_each(begin(dilation_grid_size), end(dilation_grid_size),
           [](auto& d) { d += 2; });

  auto dilation_buffer =
    vector<DilationCell>(LinearSize(dilation_grid_size), DilationCell::kBackground);
  auto dilation_grid =
    Grid<DilationCell, N>(dilation_grid_size, dilation_buffer.front());

  // Set foreground.
  for (auto const& index : indices) {
    assert(Inside(index, grid_size));
    auto dilation_index = index;
    for_each(begin(dilation_index), end(dilation_index),
             [](auto& d) { d += int32_t{1}; });
    dilation_grid.Cell(dilation_index) = DilationCell::kForeground;
  }

  // Add dilated cell indices.
  auto const dilation_neighbor_offsets = VertexNeighborOffsets<N>();
  auto dilation_indices = vector<array<int32_t, N>>();
  for (auto const& grid_index : indices) {
    assert(Inside(grid_index, grid_size));
    auto dilation_index = grid_index;
    for_each(begin(dilation_index), end(dilation_index),
             [](auto& d) { d += int32_t{1}; });
    assert(dilation_grid.Cell(dilation_index) == DilationCell::kForeground);

    for (auto const dilation_neighbor_offset : dilation_neighbor_offsets) {
      auto neighbor_index = dilation_index;
      for (auto i = size_t{0}; i < N; ++i) {
        neighbor_index[i] += dilation_neighbor_offset[i];
      }

      if (dilation_grid.Cell(neighbor_index) == DilationCell::kBackground) {
        dilation_grid.Cell(neighbor_index) = DilationCell::kDilated;
        dilation_indices.push_back(neighbor_index);
      }
    }
  }

  auto const dilation_bands_neighbor_offsets = FaceNeighborOffsets<N>();
  auto connected_dilation_components = vector<vector<array<int32_t, N>>>();
  ConnectedComponents(
    dilation_indices,
    dilation_grid_size,
    begin(dilation_bands_neighbor_offsets),
    end(dilation_bands_neighbor_offsets),
    &connected_dilation_components);

  for (auto const& dilation_component : connected_dilation_components) {
    auto dilation_band = vector<array<int32_t, N>>();
    for (auto const& dilation_index : dilation_component) {
      auto grid_index = dilation_index;
      for_each(begin(grid_index), end(grid_index),
               [](auto& d) { d -= int32_t{1}; });
      if (Inside(grid_index, grid_size)) {
        dilation_band.push_back(grid_index);
      }
    }

    if (!dilation_band.empty()) {
      dilation_bands->push_back(dilation_band);
    }
  }
}


template<std::size_t N> inline
std::array<std::pair<std::int32_t, std::int32_t>, N> BoundingBox(
  std::vector<std::array<std::int32_t, N>> const& indices)
{
  using namespace std;

  if (indices.empty()) {
    throw runtime_error("cannot compute bounding box from empty indices");
  }

  array<pair<int32_t, int32_t>, N> bbox;

  for (auto i = size_t{0}; i < N; ++i) {
    bbox[i].first = numeric_limits<int32_t>::max();
    bbox[i].second = numeric_limits<int32_t>::min();
  }

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
  auto hyper_volume = size_t{1};
  for (auto i = size_t{0}; i < N; ++i) {
    assert(bbox[i].first <= bbox[i].second);
    hyper_volume *= (bbox[i].second - bbox[i].first + 1);
  }
  return hyper_volume;
}


template<std::size_t N> inline
void InitialUnsignedNarrowBand(
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>>* narrow_band_indices)
{
  using namespace std;

  assert(!frozen_indices.empty());
  assert(narrow_band_indices != nullptr);
  narrow_band_indices->clear();

  enum class NarrowBandCell : uint8_t {
    kBackground = uint8_t{0},
    kFrozen,
    kNarrowBand
  };

  auto narrow_band_buffer =
    vector<NarrowBandCell>(LinearSize(grid_size), NarrowBandCell::kBackground);
  auto narrow_band_grid = Grid<NarrowBandCell, N>(grid_size, narrow_band_buffer);

  // Set frozen cells.
  for (auto const& frozen_index : frozen_indices) {
    assert(Inside(frozen_index, narrow_band_grid.size()));
    narrow_band_grid.Cell(frozen_index) = NarrowBandCell::kFrozen;
  }

  // Find face neighbors of frozen cells, which then become
  // the initial narrow band.
  auto const kNeighborOffsets = array<int32_t, 2>{{-1, 1}};
  for (auto const& frozen_index : frozen_indices) {
    for (auto i = size_t{0}; i < N; ++i) {
      for (auto const neighbor_offset : kNeighborOffsets) {
        auto neighbor_index = frozen_index;
        neighbor_index[i] += neighbor_offset;

        if (Inside(neighbor_index, narrow_band_grid.size())) {
          auto& neighbor_cell = narrow_band_grid.Cell(neighbor_index);
          if (neighbor_cell == NarrowBandCell::kBackground) {
            neighbor_cell = NarrowBandCell::kNarrowBand;
            narrow_band_indices->push_back(neighbor_index);
          }
        }
      }
    }
  }
}


template<std::size_t N> inline
void InitialSignedNarrowBands(
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>>* inside_narrow_band_indices,
  std::vector<std::array<std::int32_t, N>>* outside_narrow_band_indices)
{
  using namespace std;

  assert(inside_narrow_band_indices != nullptr);
  assert(outside_narrow_band_indices != nullptr);
  inside_narrow_band_indices->clear();
  outside_narrow_band_indices->clear();

  auto connected_components = vector<vector<array<int32_t, N>>>();
  auto const cc_neighbor_offsets = VertexNeighborOffsets<N>();
  ConnectedComponents(
    frozen_indices,
    grid_size,
    begin(cc_neighbor_offsets),
    end(cc_neighbor_offsets),
    &connected_components);
  assert(!connected_components.empty());

  /*
  // TODO: check for contained components!
  auto cc_bbox = vector<pair<size_t, array<pair<int32_t, int32_t>, N>>>();
  for (auto i = size_t{0}; i < connected_components.size(); ++i) {
    cc_bbox.push_back({i, BoundingBox(connected_components[i])});
  }
  sort(
    begin(cc_bbox),
    end(cc_bbox),
    [](auto const lhs, auto const rhs) {
      return HyperVolume(lhs.second) > HyperVolume(rhs.second);
    });
  for (auto i = size_t{1}; i < connected_components.size(); ++i) {
    if (Contains(cc_bbox[0].second, cc_bbox[i].second)) {
      throw runtime_error("contained connected component");
    }
  }
  */

  enum class NarrowBandCell : uint8_t {
    kBackground = uint8_t{0},
    kFrozen,
    kNarrowBand
  };

  auto narrow_band_buffer =
    vector<NarrowBandCell>(LinearSize(grid_size), NarrowBandCell::kBackground);
  auto narrow_band_grid =
    Grid<NarrowBandCell, N>(grid_size, narrow_band_buffer.front());

  // Set frozen cells.
  for (auto const& frozen_index : frozen_indices) {
    assert(Inside(frozen_index, narrow_band_grid.size()));
    narrow_band_grid.Cell(frozen_index) = NarrowBandCell::kFrozen;
  }

  for (auto const& connected_component : connected_components) {
    auto dilation_bands = vector<vector<array<int32_t, N>>>();
    DilationBands(connected_component, grid_size, &dilation_bands);
    assert(!dilation_bands.empty());
    if (dilation_bands.size() == 1) {
      throw runtime_error("open connected component");
#if 0
      if (connected_component.size() == 1) {

      }
      else {
      }
#endif
    }
    else {
      auto dilation_band_areas = vector<pair<size_t, size_t>>();
      for (auto i = size_t{0}; i < dilation_bands.size(); ++i) {
        dilation_band_areas.push_back(
          {i, HyperVolume(BoundingBox(dilation_bands[i]))});
      }
      sort(
        begin(dilation_band_areas),
        end(dilation_band_areas),
        [](auto const lhs, auto const rhs) {
          return lhs.second > rhs.second;
        });

      const auto kNeighborOffsets = array<int32_t, 2>{{-1, 1}};

      // Outer dilation bands of several connected components may overlap.
      auto const& outer_dilation_band =
        dilation_bands[dilation_band_areas[0].first];
      for (auto const& dilation_index : outer_dilation_band) {
        assert(Inside(dilation_index, narrow_band_grid.size()));
        assert(narrow_band_grid.Cell(dilation_index) != NarrowBandCell::kFrozen);
        if (narrow_band_grid.Cell(dilation_index) == NarrowBandCell::kBackground) {
          auto i = size_t{0};
          auto frozen_neighbor_found = false;
          while (i < N && !frozen_neighbor_found) {
            auto j = size_t{0};
            while (j < kNeighborOffsets.size() && !frozen_neighbor_found) {
              auto neighbor_index = dilation_index;
              neighbor_index[i] += kNeighborOffsets[j];
              if (narrow_band_grid.Cell(neighbor_index) == NarrowBandCell::kFrozen) {
                frozen_neighbor_found = true;
                narrow_band_grid.Cell(dilation_index) = NarrowBandCell::kNarrowBand;
                outside_narrow_band_indices->push_back(dilation_index);
              }
              ++j;
            }
            ++i;
          }
        }
      }

      // Inner dilation bands cannot overlap for difference connected components.
      for (auto k = size_t{1}; k < dilation_bands.size(); ++k) {
        auto const& inner_dilation_band =
          dilation_bands[dilation_band_areas[k].first];
        for (auto const& dilation_index : inner_dilation_band) {
          assert(Inside(dilation_index, narrow_band_grid.size()));
          assert(narrow_band_grid.Cell(dilation_index) == NarrowBandCell::kBackground);
          auto i = size_t{0};
          auto frozen_neighbor_found = false;
          while (i < N && !frozen_neighbor_found) {
            auto j = size_t{0};
            while (j < kNeighborOffsets.size() && !frozen_neighbor_found) {
              auto neighbor_index = dilation_index;
              neighbor_index[i] += kNeighborOffsets[j];
              if (narrow_band_grid.Cell(neighbor_index) == NarrowBandCell::kFrozen) {
                frozen_neighbor_found = true;
                narrow_band_grid.Cell(dilation_index) = NarrowBandCell::kNarrowBand;
                inside_narrow_band_indices->push_back(dilation_index);
              }
              ++j;
            }
            ++i;
          }
        }
      }
    }
  }
}


//! Set frozen cell state and distance on respective grids.
//! Throws std::runtime_error if:
//! - Frozen indices are empty, or
//! - Not the same number of indices and distances, or
//! - Any index is outside the grid (grids are asserted to have same size), or
//! - Any duplicate indices, or
//! - Any distance does not pass the predicate test, or
//! - The whole grid is frozen.
template <typename T, std::size_t N, typename D> inline
void InitializeFrozenCells(
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances,
  T const multiplier,
  D const distance_predicate,
  Grid<T, N>* const distance_grid,
  Grid<MarchingCellState, N>* const state_grid)
{
  using namespace std;

  assert(distance_grid != nullptr);
  assert(state_grid != nullptr);
  assert(distance_grid->size() == state_grid->size());

  if (frozen_indices.empty()) {
    throw runtime_error("empty frozen indices");
  }

  if (frozen_indices.size() != frozen_distances.size()) {
    throw runtime_error("frozen indices/distances size mismatch");
  }

  for (auto i = size_t{0}; i < frozen_indices.size(); ++i) {
    auto const index = frozen_indices[i];
    auto const distance = frozen_distances[i];
    if (!Inside(index, distance_grid->size())) {
      auto ss = stringstream();
      ss << "frozen index outside grid: " << ToString(index);
      throw runtime_error(ss.str());
    }

    auto& state_cell = state_grid->Cell(index);
    if (state_cell != MarchingCellState::kFar) {
      auto ss = stringstream();
      ss << "duplicate frozen index: " << ToString(index);
      throw runtime_error(ss.str());
    }
    state_cell = MarchingCellState::kFrozen;

    if (!distance_predicate(distance)) {
      auto ss = stringstream();
      ss << "invalid frozen distance: " << distance;
      throw runtime_error(ss.str());
    }
    distance_grid->Cell(index) = multiplier * distance;
  }

  // Here we know that all frozen indices are unique and inside the grid.
  if (frozen_indices.size() == LinearSize(distance_grid->size())) {
    throw std::runtime_error("whole grid frozen");
  }
}


//! Estimate distances for the initial narrow band and insert
//! indices into narrow band data structure.
template <typename T, std::size_t N, typename E> inline
void InitializeNarrowBand(
  std::vector<std::array<std::int32_t, N>> const& narrow_band_indices,
  E const& eikonal_solver,
  Grid<T, N>* const distance_grid,
  Grid<MarchingCellState, N>* const state_grid,
  NarrowBandStore<T, N>* const narrow_band)
{
  using namespace std;

  assert(distance_grid != nullptr);
  assert(state_grid != nullptr);
  assert(narrow_band != nullptr);
  assert(narrow_band->empty());

  for (auto const& narrow_band_index : narrow_band_indices) {
    assert(Inside(narrow_band_index, distance_grid->size()));
    assert(Inside(narrow_band_index, state_grid->size()));
    assert(distance_grid->Cell(narrow_band_index) == numeric_limits<T>::max());
    assert(state_grid->Cell(narrow_band_index) == MarchingCellState::kFar);

    auto const distance = eikonal_solver.Solve(
      narrow_band_index,
      *distance_grid,
      *state_grid);
    distance_grid->Cell(narrow_band_index) = distance;
    state_grid->Cell(narrow_band_index) = MarchingCellState::kNarrowBand;
    narrow_band->Insert({distance, narrow_band_index});
  }
}


template <typename T, std::size_t N, typename E> inline
void UpdateNeighbors(
  std::array<std::int32_t, N> const& index,
  E const& eikonal_solver,
  Grid<T, N>* const distance_grid,
  Grid<MarchingCellState, N>* const state_grid,
  NarrowBandStore<T, N>* const narrow_band)
{
  using namespace std;

  assert(distance_grid != nullptr);
  assert(state_grid != nullptr);
  assert(narrow_band != nullptr);
  assert(distance_grid->size() == state_grid->size());
  assert(Inside(index, distance_grid->size()));
  assert(Inside(index, state_grid->size()));
  assert(state_grid->Cell(index) == MarchingCellState::kFrozen);

  const auto kNeighborOffsets = array<int32_t, 2>{{-1, 1}};
  for (auto i = size_t{0}; i < N; ++i) {
    for (auto const neighbor_offset : kNeighborOffsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;

      if (Inside(neighbor_index, distance_grid->size())) {
        // Update the narrow band.
        auto& neighbor_state = state_grid->Cell(neighbor_index);
        switch (neighbor_state) {
        case MarchingCellState::kFar:
          {
            auto const distance = eikonal_solver.Solve(
              neighbor_index,
              *distance_grid,
              *state_grid);
            distance_grid->Cell(neighbor_index) = distance;
            neighbor_state = MarchingCellState::kNarrowBand;
            narrow_band->Insert({distance, neighbor_index});
          }
          break;
        case MarchingCellState::kNarrowBand:
          {
            auto& neighbor_distance = distance_grid->Cell(neighbor_index);
            auto const new_neighbor_distance = eikonal_solver.Solve(
              neighbor_index,
              *distance_grid,
              *state_grid);
            if (new_neighbor_distance < neighbor_distance) {
              narrow_band->DecreaseDistance(neighbor_index, new_neighbor_distance);
              neighbor_distance = new_neighbor_distance;
            }
          }
          break;
        case MarchingCellState::kFrozen:
          // If neighbor cell is frozen do nothing to it!
          break;
        }
      }
    }
  }
}


template <typename T, std::size_t N, typename E> inline
void MarchNarrowBand(
  E const& eikonal_solver,
  NarrowBandStore<T, N>* const narrow_band,
  Grid<T, N>* const distance_grid,
  Grid<MarchingCellState, N>* const state_grid)
{
  using namespace std;

  assert(distance_grid != nullptr);
  assert(state_grid != nullptr);
  assert(narrow_band != nullptr);

  while (!narrow_band->empty()) {
    // Take smallest distance from narrow band and freeze it.
    auto const narrow_band_cell = narrow_band->Pop();
    auto const distance = narrow_band_cell.first;
    auto const index = narrow_band_cell.second;

    assert(state_grid->Cell(index) == MarchingCellState::kNarrowBand);

    distance_grid->Cell(index) = distance;
    state_grid->Cell(index) = MarchingCellState::kFrozen;

    UpdateNeighbors(
      index,
      eikonal_solver,
      distance_grid,
      state_grid,
      narrow_band);
  }
}


//! Polynomial coefficients are equivalent to array index,
//! i.e. Sum(q[i] * x^i) = 0, for i in [0, 2], or simpler
//! q[0] + q[1] * x + q[2] * x^2 = 0.
//!
//! Returns the largest real root.
//!
//! Throws a runtime_error if no real roots exist.
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
  Grid<MarchingCellState, N> const& state_grid,
  T const speed,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(Inside(index, state_grid.size()));
  assert(state_grid.Cell(index) != MarchingCellState::kFrozen);
  assert(Inside(index, distance_grid.size()));

  static auto const neighbor_offsets = array<int32_t, 2>{{-1, 1}};

  auto q = array<T, 3>{{T{-1} / Squared(speed), T{0}, T{0}}};

  // Find the smallest frozen neighbor (if any) in each dimension.
  for (auto i = size_t{0}; i < N; ++i) {
    auto neighbor_min_distance = numeric_limits<T>::max();

    for (auto const neighbor_offset : neighbor_offsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;
      if (Inside(neighbor_index, distance_grid.size()) &&
          state_grid.Cell(neighbor_index) == MarchingCellState::kFrozen) {
        neighbor_min_distance =
          min(neighbor_min_distance, distance_grid.Cell(neighbor_index));
      }
    }

    // Update quadratic coefficients for the current direction.
    // If no frozen neighbor was found that dimension does not contribute
    // to the coefficients.
    if (neighbor_min_distance < numeric_limits<T>::max()) {
      auto const inv_grid_spacing_squared = T{1} / Squared(grid_spacing[i]);
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
  Grid<MarchingCellState, N> const& state_grid,
  T const speed,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  static_assert(std::is_floating_point<T>::value,
                "scalar type must be floating point");

  assert(Inside(index, state_grid.size()));
  assert(state_grid.Cell(index) != MarchingCellState::kFrozen);
  assert(Inside(index, distance_grid.size()));

  static auto const neighbor_offsets = array<int32_t, 2>{{-1, 1}};

  auto q = array<T, 3>{{T{-1} / Squared(speed), T{0}, T{0}}};

  // Find the smallest frozen neighbor (if any) in each dimension.
  for (auto i = size_t{0}; i < N; ++i) {
    auto neighbor_min_distance = numeric_limits<T>::max();
    auto neighbor_min_distance2 = numeric_limits<T>::max();

    // Check neighbors in both directions for this dimenion.
    for (auto const neighbor_offset : neighbor_offsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;
      if (Inside(neighbor_index, distance_grid.size()) &&
          state_grid.Cell(neighbor_index) == MarchingCellState::kFrozen) {
        // Neighbor one step away is frozen.
        auto const neighbor_distance = distance_grid.Cell(neighbor_index);
        if (neighbor_distance < neighbor_min_distance) {
          neighbor_min_distance = neighbor_distance;

          // Check if neighbor two steps away is frozen and has smaller
          // distance than neighbor one step away.
          auto neighbor_index2 = neighbor_index;
          neighbor_index2[i] += neighbor_offset;
          if (Inside(neighbor_index2, distance_grid.size()) &&
              state_grid.Cell(neighbor_index2) == MarchingCellState::kFrozen) {
            auto const neighbor_distance2 = distance_grid.Cell(neighbor_index2);
            if (neighbor_distance2 <= neighbor_distance) {
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
        auto const inv_grid_spacing_squared = T{1} / Squared(grid_spacing[i]);
        q[0] += Squared(neighbor_min_distance) * inv_grid_spacing_squared;
        q[1] += T{-2} * neighbor_min_distance * inv_grid_spacing_squared;
        q[2] += inv_grid_spacing_squared;
      }
    }
  }

  return SolveEikonalQuadratic(q);
}


template <typename T, std::size_t N>
class EikonalSolverBase
{
public:
  typedef T ScalarType;
  static std::size_t const kDimension = N;

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


template <typename T, std::size_t N>
class UniformSpeedEikonalSolverBase : public EikonalSolverBase<T, N>
{
public:
  UniformSpeedEikonalSolverBase(
    std::array<T, N> const& grid_spacing,
    T const speed)
    : EikonalSolverBase<T, N>(grid_spacing)
    , speed_(speed)
  {
    ThrowIfInvalidSpeed(speed_);
  }

  T speed() const
  {
    return speed_;
  }

private:
  T const speed_;
};


template <typename T, std::size_t N>
class VaryingSpeedEikonalSolverBase : public EikonalSolverBase<T, N>
{
public:
  VaryingSpeedEikonalSolverBase(
    std::array<T, N> const& grid_spacing,
    std::array<std::size_t, N> const& speed_grid_size,
    std::vector<T> const& speed_buffer)
    : EikonalSolverBase(grid_spacing)
    , speed_grid_(speed_grid_size, speed_buffer)
  {
    ThrowIfInvalidSpeedBuffer(speed_buffer);
  }

  T speed(std::array<std::int32_t, N> const& index) const
  {
    using namespace std;

    if (!Inside(index, speed_grid_.size())) {
      throw runtime_error("index outside speed grid");
    }

    return speed_grid_.Cell(index);
  }

private:
  ConstGrid<T, N> const speed_grid_;
};

} // namespace detail


//! Holds parameters related to the grid and provides methods for solving
//! the eikonal equation for a single grid cell at a time, using information
//! about both distance and state of neighboring grid cells.
template<typename T, std::size_t N>
class UniformSpeedEikonalSolver : public detail::UniformSpeedEikonalSolverBase<T, N>
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
    detail::Grid<T, N> const& distance_grid,
    detail::Grid<detail::MarchingCellState, N> const& state_grid) const
  {
    return detail::SolveEikonal(
      index,
      distance_grid,
      state_grid,
      speed(),
      grid_spacing());
  }
};


template <typename T, std::size_t N>
class HighAccuracyUniformSpeedEikonalSolver : public detail::UniformSpeedEikonalSolverBase<T, N>
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
    detail::Grid<T, N> const& distance_grid,
    detail::Grid<detail::MarchingCellState, N> const& state_grid) const
  {
    return detail::HighAccuracySolveEikonal(
      index,
      distance_grid,
      state_grid,
      speed(),
      grid_spacing());
  }
};


template <typename T, std::size_t N>
class VaryingSpeedEikonalSolver : public detail::VaryingSpeedEikonalSolverBase<T, N>
{
public:
  VaryingSpeedEikonalSolver(
    std::array<T, N> const& grid_spacing,
    std::array<std::size_t, N> const& speed_grid_size,
    std::vector<T> const& speed_buffer)
    : detail::VaryingSpeedEikonalSolverBase<T, N>(grid_spacing, speed_grid_size, speed_buffer)
  {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) and states (@a state_grid) of other cells.
  T Solve(
    std::array<std::int32_t, N> const& index,
    detail::Grid<T, N> const& distance_grid,
    detail::Grid<detail::MarchingCellState, N> const& state_grid) const
  {
    return detail::SolveEikonal(
      index,
      distance_grid,
      state_grid,
      speed(index),
      grid_spacing());
  }
};


template <typename T, std::size_t N>
class HighAccuracyVaryingSpeedEikonalSolver : public detail::VaryingSpeedEikonalSolverBase<T, N>
{
public:
  HighAccuracyVaryingSpeedEikonalSolver(
    std::array<T, N> const& grid_spacing,
    std::array<std::size_t, N> const& speed_grid_size,
    std::vector<T>& speed_buffer)
    : detail::VaryingSpeedEikonalSolverBase<T, N>(grid_spacing, speed_grid_size, speed_buffer)
  {}

  //! Returns the distance for grid cell at @a index given the current
  //! distances (@a distance_grid) and states (@a state_grid) of other cells.
  T Solve(
    std::array<std::int32_t, N> const& index,
    detail::Grid<T, N> const& distance_grid,
    detail::Grid<detail::MarchingCellState, N> const& state_grid) const
  {
    return detail::HighAccuracySolveEikonal(
      index,
      distance_grid,
      state_grid,
      speed(index),
      grid_spacing());
  }
};


//!
template<typename T, std::size_t N, typename EikonalSolverType> inline
std::vector<T> UnsignedDistance(
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances,
  EikonalSolverType const& eikonal_solver)
{
  using namespace std;
  using namespace detail;

  typedef T DistanceType;

  static_assert(is_floating_point<DistanceType>::value,
                "distance type must be floating point");
  static_assert(N > 0, "number of dimensions must be > 0");
  static_assert(N == EikonalSolverType::kDimension,
                "mismatching eikonal solver dimension");

  auto state_buffer =
    vector<MarchingCellState>(LinearSize(grid_size), MarchingCellState::kFar);
  auto distance_buffer =
    vector<DistanceType>(LinearSize(grid_size), numeric_limits<DistanceType>::max());
  auto state_grid = Grid<MarchingCellState, N>(grid_size, state_buffer);
  auto distance_grid = Grid<DistanceType, N>(grid_size, distance_buffer);
  InitializeFrozenCells(
    frozen_indices,
    frozen_distances,
    DistanceType{1}, // Multiplier
    [](auto const d) { return !isnan(d) && d >= T{0}; },
    &distance_grid,
    &state_grid);

  auto narrow_band_indices = vector<array<int32_t, N>>{};
  InitialUnsignedNarrowBand(
    frozen_indices,
    grid_size,
    &narrow_band_indices);

  auto narrow_band = NarrowBandStore<DistanceType, N>();
  InitializeNarrowBand(
    narrow_band_indices,
    eikonal_solver,
    &distance_grid,
    &state_grid,
    &narrow_band);

  MarchNarrowBand(
    eikonal_solver,
    &narrow_band,
    &distance_grid,
    &state_grid);

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
template<typename T, std::size_t N>
std::vector<T> SignedDistance(
  std::array<std::size_t, N> const& grid_size,
  std::array<T, N> const& dx,
  T const speed,
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances)
{
  using namespace std;
  using namespace detail;

  typedef T DistanceType;
  typedef EikonalSolver<DistanceType, N> EikonalSolverType;

  static_assert(is_floating_point<DistanceType>::value,
                "distance type must be floating point");
  static_assert(N > 1, "number of dimensions must be > 1");

  ThrowIfInvalidGridSize(grid_size);
  ThrowIfInvalidGridSpacing(dx);
  ThrowIfInvalidSpeed(speed);
  ThrowIfSizeNotEqual(frozen_indices, frozen_distances);
  ThrowIfEmptyIndices(frozen_indices);
  ThrowIfIndexOutsideGrid(frozen_indices, grid_size);
  ThrowIfDuplicateIndices(frozen_indices, grid_size);
  ThrowIfWholeGridFrozen(frozen_indices, grid_size);
  ThrowIfInvalidDistance(
    frozen_distances,
    [](auto const d) { return !isnan(d); });

  auto state_buffer =
    vector<MarchingCellState>(LinearSize(grid_size), MarchingCellState::kFar);
  auto distance_buffer = vector<DistanceType>(
    LinearSize(grid_size), numeric_limits<DistanceType>::max());
  auto state_grid = Grid<MarchingCellState, N>(grid_size, state_buffer.front());
  auto distance_grid = Grid<DistanceType, N>(grid_size, distance_buffer.front());
  InitializeFrozenCells(
    frozen_indices,
    frozen_distances,
    DistanceType{-1}, // Multiplier
    &distance_grid,
    &state_grid);

  auto inside_narrow_band_indices = vector<array<int32_t, N>>();
  auto outside_narrow_band_indices = vector<array<int32_t, N>>();
  InitialSignedNarrowBands(
    frozen_indices,
    grid_size,
    &inside_narrow_band_indices,
    &outside_narrow_band_indices);

  auto const eikonal_solver = EikonalSolverType(dx, speed);
  auto narrow_band = NarrowBandStore<DistanceType, N>();
  InitializeNarrowBand(
    inside_narrow_band_indices,
    eikonal_solver,
    &distance_grid,
    &state_grid,
    &narrow_band);

  MarchNarrowBand(
    eikonal_solver,
    &narrow_band,
    &distance_grid,
    &state_grid);

  // Negate all the inside distance values and flip the frozen values
  // back to their original values.
  for_each(
    begin(distance_buffer),
    end(distance_buffer),
    [](auto& d) {
      d = d < numeric_limits<DistanceType>::max() ? d * DistanceType{-1} : d;
    });

  InitializeNarrowBand(
    outside_narrow_band_indices,
    eikonal_solver,
    &distance_grid,
    &state_grid,
    &narrow_band);

  MarchNarrowBand(
    eikonal_solver,
    &narrow_band,
    &distance_grid,
    &state_grid);

  return distance_buffer;
}

} // namespace fmm
} // namespace thinks

#endif // THINKS_FMM_FASTMARCHINGMETHOD_HPP_INCLUDED
