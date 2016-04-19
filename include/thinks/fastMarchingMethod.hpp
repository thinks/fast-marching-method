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

template<typename InputIt, typename OutputIt,
         typename UnaryPredicate, typename Modifier> inline
void conditionalCopyModified(
  InputIt src_first,
  InputIt src_last,
  OutputIt dst_first,
  UnaryPredicate pred,
  Modifier mod)
{
  while (src_first != src_last) {
    if (pred(*src_first)) {
      *dst_first = mod(*src_first);
    }
    ++dst_first;
    ++src_first;
  }
}


template<typename InputIt, typename OutputIt, typename UnaryPredicate> inline
void conditionalCopy(
  InputIt src_first,
  InputIt src_last,
  OutputIt dst_first,
  UnaryPredicate pred)
{
  while (src_first != src_last) {
    if (pred(*src_first)) {
      *dst_first = *src_first;
    }
    ++dst_first;
    ++src_first;
  }
}


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
  static_assert(std::is_floating_point<T>::value, "value must be floating point");

  return T(1) / Squared(x);
}

template<typename T, std::size_t N> inline
std::array<T, N> InverseSquared(std::array<T, N> const& a)
{
  using namespace std;

  array<T, N> r;
  transform(begin(a), end(a), begin(r),
            [](T const x) { return InverseSquared(x); });
  return r;
}


template<typename T, std::size_t N> inline
T SquaredMagnitude(std::array<T, N> const& v)
{
  using namespace std;

  auto sm = T(0);
  for (auto i = size_t{0}; i < N; ++i) {
    sm += v[i] * v[i];
  }
  return sm;
}


template<std::size_t N> inline
bool Inside(
  std::array<std::int32_t, N> const& index,
  std::array<std::size_t, N> const& size)
{
  using namespace std;

  for (size_t i = 0; i < N; ++i) {
    // Cast is safe since we check that index[i] is greater than or
    // equal to zero first.
    if (!(0 <= index[i] && static_cast<size_t>(index[i]) < size[i])) {
      return false;
    }
  }
  return true;
}


template<std::size_t N> inline
void ThrowIfInvalidGridSize(std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  if (find_if(begin(grid_size), end(grid_size),
              [](auto const x) { return x == size_t{0}; }) != end(grid_size)) {
    throw runtime_error("invalid grid size");
  }
}


template<typename T, std::size_t N> inline
void ThrowIfInvalidSpacing(std::array<T, N> const& dx)
{
  using namespace std;

  if (find_if(begin(dx), end(dx),
              [](auto const x) { return x <= T(0); }) != end(dx)) {
    throw runtime_error("invalid spacing");
  }
}


template<typename T> inline
void ThrowIfInvalidSpeed(T const speed)
{
  using namespace std;

  if (speed <= T(0)) {
    auto ss = stringstream();
    ss << "invalid speed: " << speed;
    throw runtime_error(ss.str());
  }
}


template<typename U, typename V> inline
void ThrowIfSizeNotEqual(U const& u, V const& v)
{
  if (u.size() != v.size()) {
    throw std::runtime_error("size mismatch");
  }
}


template<std::size_t N> inline
void ThrowIfEmptyIndices(
  std::vector<std::array<std::int32_t, N>> const& indices)
{
  if (indices.empty()) {
    throw std::runtime_error("empty indices");
  }
}


template<std::size_t N> inline
void ThrowIfIndexOutsideGrid(
  std::vector<std::array<std::int32_t, N>> const& indices,
  std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  for (auto const& index : indices) {
    if (!Inside(index, grid_size)) {
      auto ss = stringstream();
      ss << "invalid index: [";
      for (auto i = size_t{0}; i < N; ++i) {
        ss << index[i];
        if (i != (N - 1)) {
          ss << ", ";
        }
      }
      ss << "]";
      throw runtime_error(ss.str());
    }
  }
}


template<std::size_t N> inline
void ThrowIfDuplicateIndices(
  std::vector<std::array<std::int32_t, N>> const& indices,
  std::array<std::size_t, N> const& grid_size)
{
  using namespace std;

  auto index_buffer = vector<uint8_t>(LinearSize(grid_size), 0);
  auto index_grid = Grid<uint8_t, N>(grid_size, index_buffer.front());

  for (auto const& index : indices) {
    auto& cell = index_grid.cell(index);
    if (cell == 1) {
      auto ss = stringstream();
      ss << "duplicate index: [";
      for (auto i = size_t{0}; i < N; ++i) {
        ss << index[i];
        if (i != (N - 1)) {
          ss << ", ";
        }
      }
      ss << "]";
      throw runtime_error(ss.str());
    }
    cell = 1;
  }
}


template<std::size_t N> inline
void ThrowIfWholeGridFrozen(
  std::vector<std::array<std::int32_t, N>> const& indices,
  std::array<std::size_t, N> const& grid_size)
{
  // Assumes indices are unique and inside grid.
  if (indices.size() == LinearSize(grid_size)) {
    throw std::runtime_error("whole grid frozen");
  }
}

template<typename T, typename U> inline
void ThrowIfInvalidDistance(std::vector<T> const& distances, U const unary_pred)
{
  using namespace std;

  for (auto const distance : distances) {
    if (!unary_pred(distance)) {
      auto ss = stringstream();
      ss << "invalid distance: " << distance;
      throw runtime_error(ss.str());
    }
  }
}


template<typename T, std::size_t N>
class Grid
{
public:
  typedef T cell_type;
  typedef std::array<std::size_t, N> size_type;
  typedef std::array<std::int32_t, N> index_type;

  Grid(size_type const& size, T& cells)
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

  size_type size() const
  {
    return size_;
  }

  //! Returns a reference to the cell at @a index. No range checking!
  cell_type& cell(index_type const& index)
  {
    return cells_[linearIndex_(index)];
  }

  //! Returns a const reference to the cell at @a index. No range checking!
  cell_type const& cell(index_type const& index) const
  {
    return cells_[linearIndex_(index)];
  }

private:
  //! Returns a linear (scalar) index into an array representing an
  //! N-dimensional grid for integer coordinate @a index.
  //! Note that this function does not check for integer overflow!
  std::size_t linearIndex_(index_type const& index) const
  {
    using namespace std;

    // Cast is safe since we check that index[i] is greater than or
    // equal to zero first.
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
  cell_type* const cells_;
};


template<typename T, std::size_t N>
class NarrowBandStore
{
public:
  typedef T distance_type;
  typedef std::array<std::int32_t, N> index_type;
  typedef std::pair<distance_type, index_type> value_type;

  NarrowBandStore()
  {}

  bool empty() const
  {
    assert(values_.empty() == index_to_pos_.empty());
    return values_.empty();
  }

  // O(log N)
  value_type pop()
  {
    if (empty()) {
      throw std::runtime_error("cannot pop empty narrow band store");
    }

    // Grab the top of the heap and use as return value below.
    auto const value = *values_.begin();

    // Place value from leaf level on top.
    swap_(0, values_.size() - 1);
    index_to_pos_.erase(value.second); // ~O(1), depends on hashing.
    values_.pop_back(); // O(1)
    assert(values_.size() == index_to_pos_.size());

    // Sift the new top value downwards to restore heap constraints.
    if (!empty()) {
      sift_down_(0);
    }

    return value;
  }

  void insert(value_type const& value)
  {
    if (index_to_pos_.find(value.second) != index_to_pos_.end()) {
      throw std::runtime_error("index must be unique");
    }

    // Insert value at leaf level and sift it upwards.
    auto const pos = values_.size();
    values_.push_back(value);
    index_to_pos_.insert({value.second, pos});
    sift_up_(pos);
  }

  void increase_distance(index_type const& index,
                         distance_type const new_distance)
  {
    auto const pos_iter = index_to_pos_.find(index);
    if (pos_iter == index_to_pos_.end()) {
      throw std::runtime_error("index not found");
    }

    auto& value = values_[pos_iter->second];
    if (new_distance <= value.first) {
      throw std::runtime_error("new distance must be greater than existing distance");
    }

    value.first = new_distance;
    sift_down_(pos_iter->second);
  }

  void decrease_distance(index_type const& index,
                         distance_type const new_distance)
  {
    auto const pos_iter = index_to_pos_.find(index);
    if (pos_iter == index_to_pos_.end()) {
      throw std::runtime_error("index not found");
    }

    auto& value = values_[pos_iter->second];
    if (new_distance >= value.first) {
      throw std::runtime_error("new distance must be less than existing distance");
    }

    value.first = new_distance;
    sift_up_(pos_iter->second);
  }

private:
  typedef typename std::vector<value_type>::size_type size_type_;

  void sift_up_(size_type_ const pos)
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
      swap_(pos, parent_pos);
      sift_up_(parent_pos); // Recursive!
    }
  }

  void sift_down_(size_type_ const pos)
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
      swap_(min_pos, pos);
      sift_down_(min_pos); // Recursive!
    }
  }

  static size_type_ parent_pos_(size_type_ const child_pos)
  {
    assert(child_pos > 0);
    return (child_pos - 1) / 2;
  }

  static size_type_ left_child_pos_(size_type_ const parent_pos)
  {
    return 2 * parent_pos + 1;
  }

  static size_type_ right_child_pos_(size_type_ const parent_pos)
  {
    return 2 * parent_pos + 2;
  }

  void swap_(size_type_ const pos0, size_type_ const pos1)
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
  void hashCombine_(V const& v, H const& hasher, std::size_t& seed)
  {
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  struct hash_type_
  {
    typedef index_type argument_type;
    typedef size_t result_type;

    result_type operator()(argument_type const& a) const
    {
      using namespace std;

      typedef typename argument_type::value_type value_type;
      hash<value_type> hasher;
      size_t seed = 0;
      for (auto i = size_t{0}; i < N; ++i) {
        hashCombine_(a[i], hasher, seed);
      }
      return seed;
    }
  };

  struct equal_type_
  {
    typedef bool result_type;
    typedef index_type first_argument_type;
    typedef index_type second_argument_type;

    result_type operator()(first_argument_type const& lhs,
                           second_argument_type const& rhs) const
    {
      for (auto i = size_t{0}; i < N; ++i) {
        if (lhs[i] != rhs[i]) {
          return false;
        }
      }
      return true;
    }
  };

  std::vector<value_type> values_;
  std::unordered_map<index_type, size_type_,
                     hash_type_, equal_type_> index_to_pos_;
};


enum class CellState
{
  Far = 0,
  NarrowBand,
  Frozen
};


template <typename T, std::size_t N>
class EikonalSolver
{
public:
  EikonalSolver(
    std::array<T, N> const& dx,
    T const speed)
    : inv_dx_squared_(InverseSquared(dx))
    , inv_speed_squared_(InverseSquared(speed))
  {}

  T solve(
    std::array<std::int32_t, N> const& index,
    Grid<T, N> const& distance_grid,
    Grid<CellState, N> const& state_grid) const
  {
    using namespace std;

    static const auto neighbor_offsets = array<int32_t, 2>{{-1, 1}};

    assert(state_grid.cell(index) != CellState::Frozen);
    assert(Inside(index, distance_grid.size()));

    auto q = array<T, 3>{{-inv_speed_squared_, T(0), T(0)}};
    auto n = array<T, N>{};
    fill(begin(n), end(n), numeric_limits<T>::max());

    // Find the smallest frozen neighbor (if any) in each dimension.
    for (auto i = size_t{0}; i < N; ++i) {
      for (auto const neighbor_offset : neighbor_offsets) {
        auto neighbor_index = index;
        neighbor_index[i] += neighbor_offset;

        if (Inside(neighbor_index, distance_grid.size()) &&
            state_grid.cell(neighbor_index) == CellState::Frozen) {
          n[i] = min(n[i], distance_grid.cell(neighbor_index));
        }
      }
    }
    assert(*min_element(begin(n), end(n)) < numeric_limits<T>::max());

    for (auto i = size_t{0}; i < N; ++i) {
      if (n[i] < numeric_limits<T>::max()) {
        q[0] += Squared(n[i]) * inv_dx_squared_[i];
        q[1] -= T(2) * n[i] * inv_dx_squared_[i];
        q[2] += inv_dx_squared_[i];
      }
    }

    auto const r = solveQuadratic_(q);
    assert(!isnan(r));
    assert(r >= T(0));
    assert(r >= *min_element(begin(n), end(n)));
    return r;
  }

private:
  //! Polynomial coefficients are equivalent to array index,
  //! i.e. Sum(q[i] * x^i) = 0, for i in [0, 2]
  //!
  //! Returns the largest real root of the quadratic. Assumes that
  //! real roots exist.
  template<typename T> static inline
  T solveQuadratic_(std::array<T, 3> const& q)
  {
    using namespace std;

    static_assert(is_floating_point<T>::value,
                  "quadratic coefficients must be floating point");

    assert(fabs(q[2]) > T(1e-9));

    auto const discriminant = q[1] * q[1] - T(4) * q[2] * q[0];
    if (discriminant < T(0)) {
      throw runtime_error("negative discriminant");
    }

    return (-q[1] + sqrt(discriminant)) / (T(2) * q[2]);
  }

  std::array<T, N> const inv_dx_squared_;
  T const inv_speed_squared_;
};


//! Returns base^exponent as a compile-time constant.
constexpr std::size_t pow(std::size_t const base, std::size_t const exponent)
{
  using namespace std;

  // NOTE: Cannot use loops in constexpr functions in C++11, have to use
  // recursion here.
  return exponent == 0 ? 1 : base * pow(base, exponent - 1);
}


template<std::size_t N> inline
std::array<std::array<std::int32_t, N>, pow(3, N) - 1> VertexNeighborOffsets()
{
  using namespace std;

  auto neighbor_offsets = array<array<int32_t, N>, pow(3, N) - 1>();
  auto offset = array<int32_t, N>();
  fill(begin(offset), end(offset), int32_t{-1});
  neighbor_offsets[0] = offset;
  auto offset_index = size_t{1};
  while (offset_index < neighbor_offsets.size()) {
    if (offset[N - 1] < 1) {
      ++offset[N - 1];
    }
    else {
      auto k = int32_t{N - 1};
      offset[k] = -1;
      --k;
      while (k >= 0 && offset[k] == 1) {
        offset[k] = -1;
        --k;
      }
      ++offset[k];
    }

    if (!all_of(begin(offset), end(offset), [](auto const i){ return i == 0; })) {
      neighbor_offsets[offset_index++] = offset;
    }
  }

  return neighbor_offsets;
}


template<std::size_t N> inline
std::array<std::array<std::int32_t, N>, 2 * N> FaceNeighborOffsets()
{
  using namespace std;

  array<array<int32_t, N>, 2 * N> neighbor_offsets;
  for (auto i = size_t{0}; i < N; ++i) {
    for (auto j = size_t{0}; j < N; ++j) {
      if (j == i) {
        neighbor_offsets[2 * i + 0][j] = +1;
        neighbor_offsets[2 * i + 1][j] = -1;
      }
      else {
        neighbor_offsets[2 * i + 0][j] = 0;
        neighbor_offsets[2 * i + 1][j] = 0;
      }
    }
  }
  return neighbor_offsets;
}


template<std::size_t N, typename NeighborIt> inline
void ConnectedComponents(
  std::vector<std::array<std::int32_t, N>> const& indices,
  std::array<std::size_t, N> const& grid_size,
  NeighborIt const neighbor_begin,
  NeighborIt const neighbor_end,
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
    label_grid.cell(index) = LabelCell::kForeground;
  }

  for (auto const& index : indices) {
    assert(Inside(index, label_grid.size()));
    assert(label_grid.cell(index) == LabelCell::kForeground ||
           label_grid.cell(index) == LabelCell::kLabelled);
    // Check if this index has been labelled already.
    if (label_grid.cell(index) == LabelCell::kForeground) {
      // Start a new component.
      label_grid.cell(index) = LabelCell::kLabelled;
      auto component = vector<array<int32_t, N>>();
      component.push_back(index);
      auto lifo = stack<array<int32_t, N>>();
      lifo.push(index);

      // Flood-fill current label.
      while (!lifo.empty()) {
        auto const top_index = lifo.top();
        lifo.pop();
        for (auto neighbor_iter = neighbor_begin;
             neighbor_iter != neighbor_end; ++neighbor_iter) {
          auto neighbor_index = top_index;
          for (auto i = size_t{0}; i < N; ++i) {
            neighbor_index[i] += (*neighbor_iter)[i];
          }

          if (Inside(neighbor_index, label_grid.size()) &&
              label_grid.cell(neighbor_index) == LabelCell::kForeground) {
            label_grid.cell(neighbor_index) = LabelCell::kLabelled;
            component.push_back(neighbor_index);
            lifo.push(neighbor_index);
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
    dilation_grid.cell(dilation_index) = DilationCell::kForeground;
  }

  // Add dilated cell indices.
  auto const dilation_neighbor_offsets = VertexNeighborOffsets<N>();
  auto dilation_indices = vector<array<int32_t, N>>();
  for (auto const& grid_index : indices) {
    assert(Inside(grid_index, grid_size));
    auto dilation_index = grid_index;
    for_each(begin(dilation_index), end(dilation_index),
             [](auto& d) { d += int32_t{1}; });
    assert(dilation_grid.cell(dilation_index) == DilationCell::kForeground);

    for (auto const dilation_neighbor_offset : dilation_neighbor_offsets) {
      auto neighbor_index = dilation_index;
      for (auto i = size_t{0}; i < N; ++i) {
        neighbor_index[i] += dilation_neighbor_offset[i];
      }

      if (dilation_grid.cell(neighbor_index) == DilationCell::kBackground) {
        dilation_grid.cell(neighbor_index) = DilationCell::kDilated;
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
bool Contains(
  std::array<std::pair<std::int32_t, std::int32_t>, N> const& bbox0,
  std::array<std::pair<std::int32_t, std::int32_t>, N> const& bbox1)
{
  for (auto i = size_t{0}; i < N; ++i) {
    assert(bbox0[i].first <= bbox0[i].second);
    assert(bbox1[i].first <= bbox1[i].second);
    if (bbox1[i].first < bbox0[i].first || bbox1[i].second > bbox0[i].second) {
      return false;
    }
  }
  return true;
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
  auto narrow_band_grid =
    Grid<NarrowBandCell, N>(grid_size, narrow_band_buffer.front());

  // Set frozen cells.
  for (auto const& frozen_index : frozen_indices) {
    assert(Inside(frozen_index, narrow_band_grid.size()));
    narrow_band_grid.cell(frozen_index) = NarrowBandCell::kFrozen;
  }

  // Find face neighbors of frozen cells.
  const auto kNeighborOffsets = array<int32_t, 2>{{-1, 1}};
  for (auto const& frozen_index : frozen_indices) {
    for (auto i = size_t{0}; i < N; ++i) {
      for (auto const neighbor_offset : kNeighborOffsets) {
        auto neighbor_index = frozen_index;
        neighbor_index[i] += neighbor_offset;

        if (Inside(neighbor_index, narrow_band_grid.size())) {
          auto& narrow_band_cell = narrow_band_grid.cell(neighbor_index);
          if (narrow_band_cell == NarrowBandCell::kBackground) {
            narrow_band_cell = NarrowBandCell::kNarrowBand;
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
    narrow_band_grid.cell(frozen_index) = NarrowBandCell::kFrozen;
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
        assert(narrow_band_grid.cell(dilation_index) != NarrowBandCell::kFrozen);
        if (narrow_band_grid.cell(dilation_index) == NarrowBandCell::kBackground) {
          auto i = size_t{0};
          auto frozen_neighbor_found = false;
          while (i < N && !frozen_neighbor_found) {
            auto j = size_t{0};
            while (j < kNeighborOffsets.size() && !frozen_neighbor_found) {
              auto neighbor_index = dilation_index;
              neighbor_index[i] += kNeighborOffsets[j];
              if (narrow_band_grid.cell(neighbor_index) == NarrowBandCell::kFrozen) {
                frozen_neighbor_found = true;
                narrow_band_grid.cell(dilation_index) = NarrowBandCell::kNarrowBand;
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
          assert(narrow_band_grid.cell(dilation_index) == NarrowBandCell::kBackground);
          auto i = size_t{0};
          auto frozen_neighbor_found = false;
          while (i < N && !frozen_neighbor_found) {
            auto j = size_t{0};
            while (j < kNeighborOffsets.size() && !frozen_neighbor_found) {
              auto neighbor_index = dilation_index;
              neighbor_index[i] += kNeighborOffsets[j];
              if (narrow_band_grid.cell(neighbor_index) == NarrowBandCell::kFrozen) {
                frozen_neighbor_found = true;
                narrow_band_grid.cell(dilation_index) = NarrowBandCell::kNarrowBand;
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


//! Set frozen cell state and distance.
template <typename T, std::size_t N> inline
void InitializeFrozenCells(
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances,
  T const multiplier,
  Grid<T, N>* const distance_grid,
  Grid<CellState, N>* const state_grid)
{
  using namespace std;

  assert(frozen_indices.size() == frozen_distances.size());
  assert(distance_grid != nullptr);
  assert(state_grid != nullptr);

  for (auto i = size_t{0}; i < frozen_indices.size(); ++i) {
    assert(Inside(frozen_indices[i], distance_grid->size()));
    assert(Inside(frozen_indices[i], state_grid->size()));
    assert(state_grid->cell(frozen_indices[i]) == CellState::Far);

    distance_grid->cell(frozen_indices[i]) = multiplier * frozen_distances[i];
    state_grid->cell(frozen_indices[i]) = CellState::Frozen;
  }
}


template <typename T, std::size_t N> inline
void InitializeNarrowBand(
  std::vector<std::array<std::int32_t, N>> const& narrow_band_indices,
  EikonalSolver<T, N> const& eikonal_solver,
  Grid<T, N>* const distance_grid,
  Grid<CellState, N>* const state_grid,
  NarrowBandStore<T, N>* const narrow_band)
{
  using namespace std;

  assert(narrow_band != nullptr);
  assert(narrow_band->empty());

  for (auto const& narrow_band_index : narrow_band_indices) {
    assert(Inside(narrow_band_index, distance_grid->size()));
    assert(Inside(narrow_band_index, state_grid->size()));
    assert(state_grid->cell(narrow_band_index) == CellState::Far);

    auto const distance = eikonal_solver.solve(
      narrow_band_index,
      *distance_grid,
      *state_grid);
    distance_grid->cell(narrow_band_index) = distance;
    state_grid->cell(narrow_band_index) = CellState::NarrowBand;
    narrow_band->insert({distance, narrow_band_index});
  }
}


template <typename T, std::size_t N> inline
void UpdateNeighbors(
  std::array<std::int32_t, N> const& index,
  EikonalSolver<T, N> const& eikonal_solver,
  Grid<T, N>* const distance_grid,
  Grid<CellState, N>* const state_grid,
  NarrowBandStore<T, N>* const narrow_band)
{
  using namespace std;

  static const auto neighbor_offsets = array<int32_t, 2>{{-1, 1}};

  assert(distance_grid != nullptr);
  assert(state_grid != nullptr);
  //assert same size!!
  assert(narrow_band != nullptr);
  assert(Inside(index, distance_grid->size()));
  assert(Inside(index, state_grid->size()));
  assert(state_grid->cell(index) == CellState::Frozen);

  for (auto i = size_t{0}; i < N; ++i) {
    for (auto const neighbor_offset : neighbor_offsets) {
      auto neighbor_index = index;
      neighbor_index[i] += neighbor_offset;

      if (Inside(neighbor_index, distance_grid->size())) {
        // Update the narrow band.
        auto& neighbor_state = state_grid->cell(neighbor_index);
        switch (neighbor_state) {
        case CellState::Far:
          {
            auto const distance = eikonal_solver.solve(
              neighbor_index,
              *distance_grid,
              *state_grid);
            distance_grid->cell(neighbor_index) = distance;
            neighbor_state = CellState::NarrowBand;
            narrow_band->insert(
              {
                distance,
                neighbor_index
              });
          }
          break;
        case CellState::NarrowBand:
          {
            auto& neighbor_distance = distance_grid->cell(neighbor_index);
            auto const new_neighbor_distance = eikonal_solver.solve(
              neighbor_index,
              *distance_grid,
              *state_grid);
            if (new_neighbor_distance < neighbor_distance) {
              narrow_band->decrease_distance(neighbor_index, new_neighbor_distance);
              neighbor_distance = new_neighbor_distance;
            }
          }
          break;
        // If neighbor cell is frozen do nothing to it!
        }
      }
    }
  }
}


template <typename T, std::size_t N> inline
void MarchNarrowBand(
  EikonalSolver<T, N> const& eikonal_solver,
  NarrowBandStore<T, N>* const narrow_band,
  Grid<T, N>* const distance_grid,
  Grid<CellState, N>* const state_grid)
{
  using namespace std;

  assert(distance_grid != nullptr);
  assert(state_grid != nullptr);
  assert(narrow_band != nullptr);

  while (!narrow_band->empty()) {
    // Take smallest distance from narrow band and freeze it.
    auto const narrow_band_cell = narrow_band->pop();
    auto const distance = narrow_band_cell.first;
    auto const index = narrow_band_cell.second;

    assert(state_grid->cell(index) == CellState::NarrowBand);

    distance_grid->cell(index) = distance;
    state_grid->cell(index) = CellState::Frozen;

    UpdateNeighbors(
      index,
      eikonal_solver,
      distance_grid,
      state_grid,
      narrow_band);
  }
}

} // namespace detail


//!
template<typename T, std::size_t N> inline
std::vector<T> unsignedDistance(
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
  static_assert(N > 0, "number of dimensions must be > 0");

  ThrowIfInvalidGridSize(grid_size);
  ThrowIfInvalidSpacing(dx);
  ThrowIfInvalidSpeed(speed);
  ThrowIfSizeNotEqual(frozen_indices, frozen_distances);
  ThrowIfEmptyIndices(frozen_indices);
  ThrowIfIndexOutsideGrid(frozen_indices, grid_size);
  ThrowIfDuplicateIndices(frozen_indices, grid_size);
  ThrowIfWholeGridFrozen(frozen_indices, grid_size);
  ThrowIfInvalidDistance(
    frozen_distances,
    [](auto const d) { return !isnan(d) && d >= T(0); });

  auto state_buffer = vector<CellState>(LinearSize(grid_size), CellState::Far);
  auto distance_buffer = vector<DistanceType>(
    LinearSize(grid_size), numeric_limits<DistanceType>::max());
  auto state_grid = Grid<CellState, N>(grid_size, state_buffer.front());
  auto distance_grid = Grid<DistanceType, N>(grid_size, distance_buffer.front());
  InitializeFrozenCells(
    frozen_indices,
    frozen_distances,
    DistanceType{1}, // Multiplier
    &distance_grid,
    &state_grid);

  auto narrow_band_indices = vector<array<int32_t, N>>();
  InitialUnsignedNarrowBand(
    frozen_indices,
    grid_size,
    &narrow_band_indices);

  auto const eikonal_solver = EikonalSolverType(dx, speed);
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
//!                      Euclidean distance.
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
std::vector<T> signedDistance(
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
  static_assert(N > 0, "number of dimensions must be > 0");

  ThrowIfInvalidGridSize(grid_size);
  ThrowIfInvalidSpacing(dx);
  ThrowIfInvalidSpeed(speed);
  ThrowIfSizeNotEqual(frozen_indices, frozen_distances);
  ThrowIfEmptyIndices(frozen_indices);
  ThrowIfIndexOutsideGrid(frozen_indices, grid_size);
  ThrowIfDuplicateIndices(frozen_indices, grid_size);
  ThrowIfWholeGridFrozen(frozen_indices, grid_size);
  ThrowIfInvalidDistance(
    frozen_distances,
    [](auto const d) { return !isnan(d); });

  auto state_buffer = vector<CellState>(LinearSize(grid_size), CellState::Far);
  auto distance_buffer = vector<DistanceType>(
    LinearSize(grid_size), numeric_limits<DistanceType>::max());
  auto state_grid = Grid<CellState, N>(grid_size, state_buffer.front());
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
