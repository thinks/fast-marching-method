#ifndef THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED
#define THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>


namespace thinks {
namespace fmm {
namespace detail {

//! Returns the product of the elements in array @a a.
//! Note: Not checking for overflow here!
template<std::size_t N> inline
std::size_t linearSize(std::array<std::size_t, N> const& a)
{
  using namespace std;
  return accumulate(begin(a), end(a), 1, multiplies<size_t>());
}


template<typename T> inline constexpr
T squared(T const x)
{
  return x * x;
}


template<typename T> inline constexpr
T inverseSquared(T const x)
{
  static_assert(std::is_floating_point<T>::value, "value must be floating point");
  return T(1) / squared(x);
}

template<typename T, std::size_t N> inline
std::array<T, N> inverseSquared(std::array<T, N> const& a)
{
  using namespace std;
  array<T, N> r;
  transform(begin(a), end(a), begin(r),
            [](T const x) { return inverseSquared(x); });
  return r;
}


template<typename U, typename V> inline
void throwIfSizeNotEqual(U const& u, V const& v)
{
  if (u.size() != v.size()) {
    throw std::runtime_error("size mismatch");
  }
}


template<typename U, typename V, typename W> inline
void throwIfSizeNotEqual(U const& u, V const& v, W const& w)
{
  auto const u_size = u.size();
  if (u_size != v.size() || u_size != w.size()) {
    throw std::runtime_error("size mismatch");
  }
}


template<std::size_t N>
class Index
{
public:
  typedef std::int32_t value_type;

  static constexpr std::size_t size()
  {
    return N;
  }

  explicit constexpr
  Index(std::array<value_type, N> const& v)
    : v_(v)
  {}

  //! No range checking!
  value_type& operator[](std::size_t const i)
  {
    return v_[i];
  }

  //! No range checking!
  value_type const& operator[](std::size_t const i) const
  {
    return v_[i];
  }

  std::array<value_type, N> toArray() const
  {
    return v_;
  }

private:
  std::array<value_type, N> v_;
};


template<std::size_t N> inline
Index<N> operator+(Index<N> const& lhs, Index<N> const& rhs)
{
  auto r = lhs;
  for (std::size_t i = 0; i < N; ++i) {
    r[i] += rhs[i];
  }
  return r;
}

template<> inline
Index<2> operator+(Index<2> const& lhs, Index<2> const& rhs)
{
  return Index<2>{{{lhs[0] + rhs[0], lhs[1] + rhs[1]}}};
}

template<> inline
Index<3> operator+(Index<3> const& lhs, Index<3> const& rhs)
{
  return Index<3>{{{lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]}}};
}


template<std::size_t N>
struct Neighborhood
{
  static constexpr std::array<Index<N>, 2 * N> offsets()
  {
    using namespace std;

    array<Index<N>, 2 * N> n;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        if (j == i) {
          n[2 * i + 0][j] = +1;
          n[2 * i + 1][j] = -1;
        }
        else {
          n[2 * i + 0][j] = 0;
          n[2 * i + 1][j] = 0;
        }
      }
    }
    return n;
  }
};

template<>
struct Neighborhood<2>
{
  static constexpr std::array<Index<2>, 4> offsets()
  {
    return {{
      Index<2>{{{+1,  0}}},
      Index<2>{{{-1,  0}}},
      Index<2>{{{ 0, +1}}},
      Index<2>{{{ 0, -1}}}
    }};
  }
};

template<>
struct Neighborhood<3>
{
  static constexpr std::array<Index<3>, 6> offsets()
  {
    return {{
      Index<3>{{{+1,  0,  0}}},
      Index<3>{{{-1,  0,  0}}},
      Index<3>{{{ 0, +1,  0}}},
      Index<3>{{{ 0, -1,  0}}},
      Index<3>{{{ 0,  0, +1}}},
      Index<3>{{{ 0,  0, -1}}}
    }};
  }
};


template<typename T, typename H> inline
void hashCombine(T const& v, H const& hasher, std::size_t& seed)
{
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


//! Polynomial coefficients are equivalent to array index,
//! i.e. Sum(coefficients[i] * x^i) = 0, for i in [0, 2]
//!
//! Returns the real roots of the quadratic if any exist, otherwise NaN.
//! If there are two roots the larger one is the first of the pair.
template<typename T> inline
std::pair<T, T> solveQuadratic(const std::array<T, 3> coefficients)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "quadratic coefficients must be floating point");

  T const eps = T(1e-9);

  T const c = coefficients[0];
  T const b = coefficients[1];
  T const a = coefficients[2];

  if (fabs(a) < eps) {
    if (fabs(b) < eps) {
      // c = 0, no solutions (or infinite solutions if c happens to be zero).
      return make_pair(numeric_limits<T>::quiet_NaN(),
                       numeric_limits<T>::quiet_NaN());
    }
    // bx + c = 0, one solution.
    return make_pair(-c / b, numeric_limits<T>::quiet_NaN());
  }

  if (fabs(b) < eps) {
    // ax^2 + c = 0
    T const r = std::sqrt(-c / a);
    return make_pair(r, -r);
  }

  T const discriminant_squared = b * b - 4 * a * c;
  if (discriminant_squared <= eps) {
    // Complex solution.
    return make_pair(numeric_limits<T>::quiet_NaN(),
                     numeric_limits<T>::quiet_NaN());
  }
  T const discriminant = std::sqrt(discriminant_squared);

  T const r0 = (b < T(0)) ?
    (-b + discriminant) / (2 * a) : // b < 0
    (-b - discriminant) / (2 * a);  // b > 0
  T const r1 = c / (a * r0);
  return make_pair(max(r0, r1), min(r0, r1));
}


template<std::size_t N> inline
bool inside(std::array<std::size_t, N> const& size, Index<N> const& index)
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
    std::size_t stride = 1;
    for (std::size_t i = 1; i < N; ++i) {
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
    size_t k = index[0];
    for (size_t i = 1; i < N; ++i) {
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
  typedef Index<N> index_type;
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

  std::vector<value_type> values_;
  std::unordered_map<index_type, size_type_> index_to_pos_;
};


enum class State
{
  Far = 0,
  NarrowBand,
  Frozen
};


template<typename T>
struct Cell
{
  T distance;
  State state;
};


template <typename T, std::size_t N>
class EikonalSolver
{
public:
  EikonalSolver(
    std::array<T, N> const& dx,
    T const speed,
    T const multiplier = T(1))
    : neighbor_offsets_(Neighborhood<N>::offsets())
    , inv_dx_squared_(inverseSquared(dx))
    , inv_speed_squared_(inverseSquared(speed))
    , multiplier_(multiplier)
  {}

  T solve(Index<N> const& index, Grid<Cell<T>, N> const& grid) const
  {
    using namespace std;

    auto q = array<T, 3>{{-inv_speed_squared_, T(0), T(0)}};

    for (size_t i = 0; i < N; ++i) {
      auto min_frozen_neighbor_distance = numeric_limits<T>::max();
      for (size_t j = 0; j < 2; ++j) {
        auto const neighbor_index = index + neighbor_offsets_[2 * i + j];
        if (inside(grid.size(), neighbor_index)) {
          auto const neighbor_cell = grid.cell(neighbor_index.toArray());
          if (neighbor_cell.state == State::Frozen) {
            assert(fabs(neighbor_cell.distance) <= fabs(grid.cell(index.toArray()).distance));
            min_frozen_neighbor_distance = min(
              min_frozen_neighbor_distance, multiplier_ * neighbor_cell.distance);
          }
        }
      }

      if (min_frozen_neighbor_distance < numeric_limits<T>::max()) {
        q[0] += squared(min_frozen_neighbor_distance) * inv_dx_squared_[i];
        q[1] += T(-2) * min_frozen_neighbor_distance * inv_dx_squared_[i];
        q[2] += inv_dx_squared_[i];
      }
    }

    auto const r = solveQuadratic(q);
    assert(!isnan(r.first));
    assert(r.first >= T(0));
    return r.first;
  }

private:
  std::array<Index<N>, 2 * N> const neighbor_offsets_;
  std::array<T, N> const inv_dx_squared_;
  T const inv_speed_squared_;
  T const multiplier_;
};


template <typename T, typename P, std::size_t N> inline
void updateNeighbors(
  EikonalSolver<T, N> const& eikonal_solver,
  Index<N> const& frozen_index,
  std::array<T, N> const& normal,
  P const pred,
  Grid<Cell<T>, N>* grid,
  NarrowBandStore<T, N>* narrow_band)
{
  assert(grid != nullptr);
  assert(narrow_band != nullptr);
  assert(grid->cell(frozen_index.toArray()).state == State::Frozen);

  // Update neighbors of frozen cell.
  for (auto const& neighbor_offset : Neighborhood<N>::offsets()) {
    auto const neighbor_index = frozen_index + neighbor_offset;
    if (inside(grid->size(), neighbor_index) && pred(normal, neighbor_offset)) {
      // Update narrow band.
      auto& neighbor_cell = grid->cell(neighbor_index.toArray());
      switch (neighbor_cell.state) {
      case State::Far:
        neighbor_cell.distance = eikonal_solver.solve(neighbor_index, *grid);
        neighbor_cell.state = State::NarrowBand;
        narrow_band->insert(
          {
            fabs(neighbor_cell.distance),
            neighbor_index
          });
        break;
      case State::NarrowBand:
        auto const new_distance = eikonal_solver.solve(neighbor_index, *grid);
        if (fabs(new_distance) < fabs(neighbor_cell.distance)) {
          narrow_band->decrease_distance(neighbor_index, fabs(new_distance));
        }
        else if (fabs(new_distance) > fabs(neighbor_cell.distance)) {
          // Could happen because of numerical instability.
          narrow_band->increase_distance(neighbor_index, new_distance);
        }
        neighbor_cell.distance = new_distance;
        break;
      // If neighbor is frozen do nothing!
      }
    }
  }
}


//! Set frozen cell state and distance.
template <typename T, typename P, std::size_t N> inline
void initializeFrozenCells(
  std::vector<std::array<typename Index<N>::value_type, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances,
  P const pred,
  Grid<Cell<T>, N>* grid)
{
  using namespace std;

  assert(grid != nullptr);

  throwIfSizeNotEqual(frozen_indices, frozen_distances);

  for (size_t i = 0; i < frozen_indices.size(); ++i) {
    auto const frozen_index = frozen_indices[i];
    if (!inside(grid->size(), Index<N>(frozen_index))) {
      throw runtime_error("frozen cell index outside grid");
    }

    auto const frozen_distance = frozen_distances[i];
    if (!pred(frozen_distance)) {
      throw runtime_error("frozen cell distance predicate failed");
    }

    auto& cell = grid->cell(frozen_index);
    assert(cell.state == State::Far);
    cell.state = State::Frozen;
    cell.distance = frozen_distance;
  }
}


template <typename T, std::size_t N, typename P> inline
std::unique_ptr<NarrowBandStore<T, N>> initializeNarrowBand(
  EikonalSolver<T, N> const& eikonal_solver,
  std::vector<std::array<typename Index<N>::value_type, N>> const& frozen_indices,
  std::vector<std::array<T, N>> const*const normals,
  P const& pred,
  Grid<Cell<T>, N>* grid)
{
  using namespace std;

  assert(grid != nullptr);

  if (normals != nullptr) {
    throwIfSizeNotEqual(*normals, frozen_indices);
  }

  auto narrow_band = make_unique<NarrowBandStore<T, N>>();

  array<T, N> dummy_normal;
  fill(begin(dummy_normal), end(dummy_normal), T(0));

  // Initialize the narrow band cells.
  for (size_t i = 0; i < frozen_indices.size(); ++i) {
    assert(inside(grid->size(), Index<N>(frozen_indices[i])));
    updateNeighbors(
      eikonal_solver,
      Index<N>(frozen_indices[i]),
      normals != nullptr ? (*normals)[i] : dummy_normal,
      pred,
      grid,
      narrow_band.get());
  }

  if (narrow_band->empty()) {
    throw runtime_error("narrow band cannot be empty after initialization");
  }

  return move(narrow_band);
}


template <typename T, std::size_t N> inline
void marchNarrowBand(
  EikonalSolver<T, N> const& eikonal_solver,
  Grid<Cell<T>, N>* grid,
  NarrowBandStore<T, N>* narrow_band,
  T const multiplier = T(1))
{
  using namespace std;

  assert(grid != nullptr);
  assert(narrow_band != nullptr);

  array<T, N> dummy_normal;
  fill(begin(dummy_normal), end(dummy_normal), T(0));

  while (!narrow_band->empty()) {
    // Take smallest distance from narrow band and freeze it.
    auto const value = narrow_band->pop();
    auto const frozen_distance = value.first;
    auto const frozen_index = value.second;
    auto& cell = grid->cell(frozen_index.toArray());
    assert(cell.state == State::NarrowBand);
    cell.distance = multiplier * frozen_distance;
    cell.state = State::Frozen;
    updateNeighbors(
      eikonal_solver,
      frozen_index,
      dummy_normal,
      [](auto const&, auto const&) {
        return true; // Always update all non-frozen neighbors while marching.
      },
      grid,
      narrow_band);
  }
}

} // namespace detail


template<typename T, std::size_t N> inline
std::vector<T> unsignedDistance(
  std::array<std::size_t, N> const& size,
  std::array<T, N> const& dx,
  T const speed,
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances)
{
  using namespace std;
  using namespace detail;

  typedef T DistanceType;
  typedef Cell<DistanceType> CellType;
  typedef Grid<CellType, N> GridType;
  typedef EikonalSolver<DistanceType, N> EikonalSolverType;

  auto cell_buffer = vector<CellType>(linearSize(size),
    CellType{numeric_limits<DistanceType>::max(), State::Far});
  auto grid = GridType(size, cell_buffer.front());

  initializeFrozenCells(
    frozen_indices,
    frozen_distances,
    [](auto const d) {
      return d >= DistanceType(0); // Input distances must be positive.
    },
    &grid);

  auto const eikonal_solver = EikonalSolverType(dx, speed);
  auto narrow_band = initializeNarrowBand<DistanceType, N>(
    eikonal_solver,
    frozen_indices,
    nullptr, // Normals.
    [](auto const&, auto const&) { return true; },
    &grid);
  marchNarrowBand(eikonal_solver, &grid, narrow_band.get());

  auto distance_buffer = vector<DistanceType>(cell_buffer.size());
  transform(begin(cell_buffer), end(cell_buffer), begin(distance_buffer),
            [](auto const& cell) { return cell.distance; });

  return distance_buffer;
}


template<typename T, std::size_t N>
std::vector<T> signedDistance(
  std::array<std::size_t, N> const& size,
  std::array<T, N> const& dx,
  T const speed,
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances,
  std::vector<std::array<T, N>> const& normals)
{
  using namespace std;
  using namespace detail;

  typedef T DistanceType;
  typedef Cell<DistanceType> CellType;
  typedef Grid<CellType, N> GridType;
  typedef EikonalSolver<DistanceType, N> EikonalSolverType;

  auto cell_buffer = vector<CellType>(linearSize(size),
    CellType{numeric_limits<DistanceType>::max(), State::Far});
  auto grid = GridType(size, cell_buffer.front());

  initializeFrozenCells(
    frozen_indices,
    frozen_distances,
    [](DistanceType const) { return true; },
    &grid);

  auto outside_eikonal_solver = EikonalSolverType(dx, speed);
  auto outside_narrow_band = initializeNarrowBand<DistanceType, N>(
    outside_eikonal_solver,
    frozen_indices,
    &normals,
    [](auto const& normal, auto const& neighbor_offset) {
      auto sum = DistanceType(0);
      for (size_t i = 0; i < N; ++i) {
        sum += normal[i] * neighbor_offset[i];
      }
      return sum >= DistanceType(0);
    },
    &grid);
  marchNarrowBand(outside_eikonal_solver, &grid, outside_narrow_band.get());

  auto inside_eikonal_solver = EikonalSolverType(dx, speed, T(-1));
  auto inside_narrow_band = initializeNarrowBand<DistanceType, N>(
    inside_eikonal_solver,
    frozen_indices,
    &normals,
    [](auto const& normal, auto const& neighbor_offset) {
      auto sum = DistanceType(0);
      for (size_t i = 0; i < N; ++i) {
        // Flip normal.
        sum += (DistanceType(-1) * normal[i]) * neighbor_offset[i];
      }
      return sum >= DistanceType(0);
    },
    &grid);
  marchNarrowBand(inside_eikonal_solver, &grid, inside_narrow_band.get(), T(-1));

  auto distance_buffer = vector<DistanceType>(cell_buffer.size());
  transform(begin(cell_buffer), end(cell_buffer), begin(distance_buffer),
            [](auto const& cell){ return cell.distance; });

  return distance_buffer;
}

} // namespace fmm
} // namespace thinks


namespace std {

// std::hash specializations for thinks::fmm::detail::Index types.

template<size_t N>
struct hash<thinks::fmm::detail::Index<N>>
{
  typedef thinks::fmm::detail::Index<N> argument_type;
  typedef size_t result_type;

  result_type operator()(argument_type const& a) const
  {
    typedef typename argument_type::value_type value_type;
    hash<value_type> hasher;
    size_t seed = 0;
    for (size_t i = 0; i < N; ++i) {
      thinks::fmm::detail::hashCombine(a[i], hasher, seed);
    }
    return seed;
  }
};

template<>
struct hash<thinks::fmm::detail::Index<2>>
{
  typedef thinks::fmm::detail::Index<2> argument_type;
  typedef size_t result_type;

  result_type operator()(argument_type const& a) const
  {
    typedef argument_type::value_type value_type;
    hash<value_type> hasher;
    size_t seed = 0;
    thinks::fmm::detail::hashCombine(a[0], hasher, seed);
    thinks::fmm::detail::hashCombine(a[1], hasher, seed);
    return seed;
  }
};

template<>
struct hash<thinks::fmm::detail::Index<3>>
{
  typedef thinks::fmm::detail::Index<3> argument_type;
  typedef size_t result_type;

  result_type operator()(argument_type const& a) const
  {
    typedef argument_type::value_type value_type;
    hash<value_type> hasher;
    size_t seed = 0;
    thinks::fmm::detail::hashCombine(a[0], hasher, seed);
    thinks::fmm::detail::hashCombine(a[1], hasher, seed);
    thinks::fmm::detail::hashCombine(a[2], hasher, seed);
    return seed;
  }
};


// std::equal_to specializations for thinks::detail::Index types.

template<size_t N>
struct equal_to<thinks::fmm::detail::Index<N>>
{
  typedef bool result_type;
  typedef thinks::fmm::detail::Index<N> first_argument_type;
  typedef thinks::fmm::detail::Index<N> second_argument_type;

  result_type operator()(first_argument_type const& lhs,
                         second_argument_type const& rhs) const
  {
    for (size_t i = 0; i < N; ++i) {
      if (lhs[i] != rhs[i]) {
        return false;
      }
    }
    return true;
  }
};

template<>
struct equal_to<thinks::fmm::detail::Index<2>>
{
  typedef bool result_type;
  typedef thinks::fmm::detail::Index<2> first_argument_type;
  typedef thinks::fmm::detail::Index<2> second_argument_type;

  result_type operator()(first_argument_type const& lhs,
                         second_argument_type const& rhs) const
  {
    return lhs[0] == rhs[0] && lhs[1] == rhs[1];
  }
};

template<>
struct equal_to<thinks::fmm::detail::Index<3>>
{
  typedef bool result_type;
  typedef thinks::fmm::detail::Index<3> first_argument_type;
  typedef thinks::fmm::detail::Index<3> second_argument_type;

  result_type operator()(first_argument_type const& lhs,
                         second_argument_type const& rhs) const
  {
    return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
  }
};

} // namespace std

#endif // THINKS_FASTMARCHINGMETHOD_HPP_INCLUDED
