#ifndef THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED
#define THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>


namespace thinks {
namespace detail {

//! Returns the product of the elements in array @a a.
template<std::size_t N> inline
std::size_t linear_size(std::array<std::size_t, N> const& a)
{
  return std::accumulate(std::begin(a), std::end(a),
                         1, std::multiplies<std::size_t>());
}


//! Generic, not implemented.
template<std::size_t D>
struct Neighborhood;

template<>
struct Neighborhood<2>
{
  static constexpr std::array<std::array<std::int32_t, 2>, 4> offsets()
  {
    return {{
        {{+1,  0}},
        {{-1,  0}},
        {{ 0, +1}},
        {{ 0, -1}}
      }};
  }
};

template<>
struct Neighborhood<3>
{
  static constexpr std::array<std::array<std::int32_t, 3>, 6> offsets()
  {
    return {{
      {{+1,  0,  0}},
      {{-1,  0,  0}},
      {{ 0, +1,  0}},
      {{ 0, -1,  0}},
      {{ 0,  0, +1}},
      {{ 0,  0, -1}}
    }};
  }
};


template<typename T, typename H> inline
void hash_combine(std::size_t& seed, T const& v, H const& hasher)
{
  seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template<typename T, std::size_t N>
struct array_hash
{
  std::size_t operator()(std::array<T, N> const& a)const
  {
    std::hash<T> hasher;
    std::size_t seed = 0;
    for (auto const& element : a) {
      hash_combine(seed, element, hasher);
    }
    return seed;
  }
};

template<typename T>
struct array_hash<T, 2>
{
  std::size_t operator()(std::array<T, 2> const& a)const
  {
    std::hash<T> hasher;
    std::size_t seed = 0;
    hash_combine(seed, a[0], hasher);
    hash_combine(seed, a[1], hasher);
    return seed;
  }
};

template<typename T>
struct array_hash<T, 3>
{
  std::size_t operator()(std::array<T, 3> const& a)const
  {
    std::hash<T> hasher;
    std::size_t seed = 0;
    hash_combine(seed, a[0], hasher);
    hash_combine(seed, a[1], hasher);
    hash_combine(seed, a[2], hasher);
    return seed;
  }
};


template<typename T, std::size_t N>
struct array_equal
{
  bool operator()(std::array<T, N> const& a, std::array<T, N> const& b) const
  {
    typedef typename std::array<T, N> size_type;
    for (size_type i = 0; i < N; ++i) {
      if (a[i] != b[i]) {
        return false;
      }
    }
    return true;
  }
};

template<typename T>
struct array_equal<T, 2>
{
  bool operator()(std::array<T, 2> const& a, std::array<T, 2> const& b) const
  {
    return a[0] == b[0] && a[1] == b[1];
  }
};

template<typename T>
struct array_equal<T, 3>
{
  bool operator()(std::array<T, 3> const& a, std::array<T, 3> const& b) const
  {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
  }
};


//! Polynomial coefficients are equivalent with array index,
//! i.e. Sum(coefficients[i] * x^i) = 0, for i in [0, 2]
//!
//! Returns the real roots of the quadratic if any exist, otherwise NaN.
template<typename T> inline
std::pair<T, T> solveQuadratic(const std::array<T, 3> coefficients)
{
  static_assert(std::is_floating_point<T>::value,
                "quadratic coefficients must be floating point");

  T const c = coefficients[0];
  T const b = coefficients[1];
  T const a = coefficients[2];

  if (std::fabs(a) < 1e-9) {
    if (std::fabs(b) < 1e-9) {
      // c = 0, no solutions (or infinite solutions if c happens to be zero).
      return std::make_pair(std::numeric_limits<T>::quiet_NaN(),
                            std::numeric_limits<T>::quiet_NaN());
    }
    // bx + c = 0, one solution.
    return std::make_pair(-c / b,
                          std::numeric_limits<T>::quiet_NaN());
  }
  // else: ax^2 + bx + c = 0

  if (std::fabs(b) < 1e-9) {
    // ax^2 + c = 0
    T const r = std::sqrt(c / a);
    return std::make_pair(r, -r);
  }

  T const discriminant_squared = b * b - 4 * a * c;
  if (discriminant_squared <= 1e-9) {
    // Complex solution.
    return std::make_pair(std::numeric_limits<T>::quiet_NaN(),
                          std::numeric_limits<T>::quiet_NaN());
  }
  T const discriminant = std::sqrt(discriminant_squared);

  T const r0 = (b < T(0)) ?
    (-b + discriminant) / (2 * a) : // b < 0
    (-b - discriminant) / (2 * a);  // b > 0
  T const r1 = c / (a * r0);
  return std::make_pair(std::max(r0, r1), std::min(r0, r1));
}


template<typename T, std::size_t N>
class Grid
{
public:
  typedef T cell_type;
  typedef std::array<std::size_t, N> index_type;

  Grid(std::array<std::size_t, N> const& size, T& cells)
    : size_(size)
    , cells_(&cells)
  {
    std::size_t stride = 1;
    for (std::size_t i = 1; i < N; ++i) {
      stride *= size_[i - 1];
      strides_[i - 1] = stride;
    }
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
  std::size_t linearIndex_(index_type const& index)
  {
    assert(size_[0] <= index[0] && index[0] < size_[0]);
    std::size_t k = index[0];
    for (std::size_t i = 1; i < N; ++i) {
      assert(size_[i] <= index[i] && index[i] < size_[i]);
      k += index[i] * strides_[i - 1];
    }
    return k;
  }

  std::array<std::size_t, N> size_;
  std::array<std::size_t, N - 1> strides_;
  cell_type * const cells_;
};


template<typename T, std::size_t N>
class NarrowBandStore
{
public:
  typedef T distance_type;
  typedef std::array<std::size_t, N> index_type;
  typedef std::pair<distance_type, index_type> value_type;
  typedef typename std::vector<value_type>::size_type size_type;

  NarrowBandStore()
  {}

#if 0
  size_type size() const
  {
    return values_.size();
  }
#endif

  bool empty() const
  {
    return values_.empty();
  }

  value_type pop()
  {
    if (empty()) {
      throw std::runtime_error("cannot pop empty narrow band store");
    }

    auto const value = *values_.begin();

    *values_.begin() = *values_.rbegin();
    values_.pop_back(); // O(1)

    if (!empty()) {
      bubble_down_(0);
    }

    return value;
  }

  void insert(value_type const& value)
  {
    values_.push_back(value);
    bubble_up_(values_.size() - 1);
  }

  //! Valid until next call to pop/insert.
  value_type const* find_index(index_type const& index) const
  {
    auto const value_iter = index_to_value_.find(index);
    return value_iter != index_to_value_.end() ?
      &values_[value_iter->second] : nullptr;
  }

  void decrease_distance(value_type const* const pos, distance_type const distance)
  {

  }

private:
#if 0
  void heapify_()
  {
    for (std::int32_t i = values_.size() / 2 - 1; i >= 0; --i) {
      bubble_down_(static_cast<std::size_t>(i));
    }
  }
#endif

  void bubble_up_(size_type const index)
  {
    if (index == 0) {
      return;
    }

    size_type const parent_index = (index - 1) / 2;
    assert(parent_index < values_.size());
    auto& parent_value = values_[parent_index];
    auto& index_value = values_[index];
    if (parent_value > index_value) {
      std::swap(parent_value, index_value);
      bubble_up_(parent_index);
    }
  }

  void bubble_down_(size_type const index)
  {
    auto const left_index = 2 * index + 1;
    auto const right_index = 2 * index + 2;
    auto const size = values_.size();

    if (left_index >= size) {
      return; // Index is a leaf.
    }

    auto const value = values_[index];
    auto const left_value = values_[left_index];
    auto min_index = index;

    if (value > left_value) {
      min_index = left_index;
    }

    if (right_index < size) {
      auto const min_value = values_[min_index];
      auto const right_value = values_[right_index];
      if (min_value > right_value) {
        min_index = right_index;
      }
    }

    if (min_index != index) {
      std::swap(values_[index], values_[min_index]);
      bubble_down_(min_index);
    }
  }

  typedef array_hash<typename index_type::value_type, N> hash_type_;
  typedef array_equal<typename index_type::value_type, N> equal_type_;

  std::vector<value_type> values_;
  std::unordered_map<index_type, size_type,
                     hash_type_, equal_type_> index_to_value_;
};

enum class State
{
  Far,
  NarrowBand,
  Frozen
};

template<typename T, std::size_t N> inline
T compute_distance(
  Grid<T, N>& distance_grid,
  std::array<std::size_t, N> const& index)
{
  // TODO: implement!
  return T(0);
}


template <typename T, std::size_t N> inline
void initializeNarrowBand(
  std::vector<std::array<std::size_t, N>> const& frozen_cell_indices,
  Grid<T, N>& distance_grid,
  Grid<State, N>& state_grid,
  NarrowBandStore<T, N>& narrow_band)
{
  if (frozen_cell_indices.empty()) {
    throw std::runtime_error("must have at least one frozen cell in order to "
                             "initialize narrow band");
  }

  for (auto const& frozen_cell_index : frozen_cell_indices) {
    state_grid.cell(frozen_cell_index) = State::Frozen;
    for (auto const& neighbor_offset : Neighborhood<N>::offsets()) {
      auto neighbor_cell_index = frozen_cell_index;
      std::size_t d = 0;
      for (auto const& offset : neighbor_offset) {
        neighbor_cell_index[d++] += offset;
      }

      // TODO: check if inside grid!!

      // Compute distance at neighbor.
      auto const distance = compute_distance(distance_grid, neighbor_cell_index);

      // Update narrow band.
      auto& neighbor_state = state_grid.cell(neighbor_cell_index);
      switch (neighbor_state) {
      case State::Far:
        neighbor_state = State::NarrowBand;
        narrow_band.insert(std::make_pair(distance, neighbor_cell_index));
        break;
      case State::NarrowBand:
        {
          auto const iter = narrow_band.find_index(neighbor_cell_index);
          assert(iter != nullptr);
          narrow_band.decrease_distance(iter, distance);
        }
        break;
      }
    }
  }
}

template <typename T, std::size_t N> inline
void marchNarrowBand(Grid<T, N>& distance_grid,
                     Grid<State, N>& state_grid,
                     NarrowBandStore<T, N>& narrow_band)
{
  while (!narrow_band.empty()) {
    auto const value = narrow_band.pop();
    auto const cell_index = value.second;
    assert(state_grid.cell(cell_index) == State::NarrowBand);
    state_grid.cell(cell_index) = State::Frozen;
    for (auto const& neighbor_offset : Neighborhood<N>::offsets()) {
      auto neighbor_cell_index = cell_index;
      std::size_t d = 0;
      for (auto const& offset : neighbor_offset) {
        neighbor_cell_index[d++] += offset;
      }

      // TODO: check if inside grid!!

      // Compute distance at neighbor.
      auto const distance = compute_distance(distance_grid, neighbor_cell_index);

      // Update narrow band.
      auto& neighbor_state = state_grid.cell(neighbor_cell_index);

      switch (neighbor_state) {
      case State::Far:
        neighbor_state = State::NarrowBand;
        narrow_band.insert(std::make_pair(distance, neighbor_cell_index));
        break;
      case State::NarrowBand:
        {
          auto const iter = narrow_band.find_index(neighbor_cell_index);
          assert(iter != nullptr);
          narrow_band.decrease_distance(iter, distance);
        }
        break;
      }
    }
  }
}

} // namespace detail

template<typename T, std::size_t N>
void fastMarchingMethod(
  std::vector<std::array<std::size_t, N>> const& frozen_cells,
  std::array<std::size_t, N> const& size,
  T* distance_buffer)
{
#if 1
  using namespace std;
  using namespace detail;

  typedef T DistanceType;
  typedef array<size_t, N> IndexType;

  assert(distance_buffer != nullptr);

  vector<State> state(linear_size(size), State::Far);
  Grid<State, N> state_grid(size, state.front());

  Grid<DistanceType, N> distance_grid(size, distance_buffer[0]);
  NarrowBandStore<DistanceType, N> narrow_band;

  initializeNarrowBand(frozen_cells, distance_grid, state_grid, narrow_band);
  marchNarrowBand(distance_grid, state_grid, narrow_band);
#endif

#if 0
  std::vector<float> v;
  v.push_back(23.f);
  v.push_back(0.f);
  v.push_back(-23.f);
  v.push_back(67.f);
  v.push_back(12.f);
  v.push_back(35.f);
  v.push_back(43.f);
  v.push_back(27.f);
  v.push_back(123.f);

  {
    detail::min_heap<float> h(v);
    while (!h.empty()) {
      std::cerr << h.pop() << std::endl;
    }
  }

  std::cerr << "-----------" << std::endl;

  {
    detail::min_heap<float> h(v);
    h.insert(-512.f);
    h.insert(512.f);

    while (!h.empty()) {
      std::cerr << h.pop() << std::endl;
    }
  }
#endif
}

} // namespace thinks

#endif // THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED
