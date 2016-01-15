#ifndef THINKS_TESTFASTMARCHINGMETHOD_HPP_INCLUDED
#define THINKS_TESTFASTMARCHINGMETHOD_HPP_INCLUDED

#include <array>
#include <cassert>
#include <utility>
#include <vector>

#include <thinks/fastMarchingMethod.hpp>


namespace thinks {
namespace fmm {
namespace test {

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
    assert(0 <= index[0] &&
      index[0] < static_cast<index_type::value_type>(size_[0]));
    std::size_t k = index[0];
    for (std::size_t i = 1; i < N; ++i) {
      assert(0 <= index[i] &&
        index[i] < static_cast<index_type::value_type>(size_[i]));
      k += index[i] * strides_[i - 1];
    }
    return k;
  }

  std::array<std::size_t, N> const size_;
  std::array<std::size_t, N - 1> strides_;
  cell_type* const cells_;
};

namespace detail {

//! Returns the product of the elements in array @a a.
//! Note: Not checking for overflow here!
template<std::size_t N> inline
std::size_t linearSize(std::array<std::size_t, N> const& a)
{
  using namespace std;
  return accumulate(begin(a), end(a), 1, multiplies<size_t>());
}


template<typename T, typename R, typename U> inline
std::vector<R> transformedVector(
  std::vector<T> const& v,
  U unary_op)
{
  using namespace std;

  auto r = vector<R>(v.size());
  transform(begin(v), end(v), begin(r), unary_op);
  return r;
}


template<std::size_t N, typename T> inline
std::array<T, N> filledArray(T const v)
{
  using namespace std;

  auto r = array<T, N>{};
  fill(begin(r), end(r), v);
  return r;
}


template<typename T, std::size_t N> inline
std::array<T, N> sub(std::array<T, N> const& u, std::array<T, N> const& v)
{
  using namespace std;

  auto r = std::array<T, N>{};
  transform(begin(u), end(u), begin(v), begin(r),
            [](auto const ux, auto const vx){ return ux - vx; });
  return r;
}


template<std::size_t N>
struct Neighborhood
{
  static std::array<std::array<std::int32_t, N>, 2 * N> offsets()
  {
    std::array<std::array<std::int32_t, N>, 2 * N> n;
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < N; ++j) {
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


template<typename T, std::size_t N> inline
std::array<T, N> gradient(
  Grid<T, N> const& grid,
  typename Grid<T, N>::index_type const& index,
  std::array<T, N> const& dx,
  std::array<std::array<std::int32_t, N>, 2 * N> const& neighbor_offsets)
{
  using namespace std;

  typedef typename Grid<T, N>::index_type IndexType;

  array<T, N> grad;
  auto const size = grid.size();

  for (size_t i = 0; i < N; ++i) {
    auto const neighbor_pos_offset = neighbor_offsets[2 * i  + 0];
    auto const neighbor_neg_offset = neighbor_offsets[2 * i  + 1];
    auto pos_index = index;
    auto neg_index = index;
    for (size_t j = 0; j < N; ++j) {
      pos_index[j] += neighbor_pos_offset[j];
      neg_index[j] += neighbor_neg_offset[j];
    }

    auto min_value = numeric_limits<T>::max();
    if (0 <= pos_index[i] && static_cast<size_t>(pos_index[i]) < size[i]) {
      min_value = grid.cell(pos_index);
    }
    if (0 <= neg_index[i] && static_cast<size_t>(neg_index[i]) < size[i]) {
      min_value = min(min_value, grid.cell(neg_index));
    }
    grad[i] = (grid.cell(index) - min_value) / dx[i];
  }
  return grad;
}


template<typename T, std::size_t N> inline
T magnitude(std::array<T, N> const& v)
{
  auto mag_squared = T{0};
  for (std::size_t i = 0; i < N; ++i) {
    mag_squared += v[i] * v[i];
  }
  return std::sqrt(mag_squared);
}


template<typename T, std::size_t N> inline
T distance(std::array<T, N> const& u, std::array<T, N> const& v)
{
  auto distance_squared = T{0};
  for (std::size_t i = 0; i < N; ++i) {
    auto const delta = u[i] - v[i];
    distance_squared += delta * delta;
  }
  return std::sqrt(distance_squared);
}


template<typename T, std::size_t N> inline
std::array<T, N> normalized(std::array<T, N> const& v)
{
  auto n = v;
  auto const mag = magnitude(v);
  for (std::size_t i = 0; i < N; ++i) {
    n[i] /= mag;
  }
  return n;
}


struct Stats
{
  double min;
  double max;
  double avg;
  double std_dev;
};


template<typename T> inline
Stats stats(std::vector<T> const& v)
{
  using namespace std;

  static_assert(is_floating_point<T>::value, "value type must be floating point");

  if (v.empty()) {
    return Stats{
      numeric_limits<T>::quiet_NaN(),
      numeric_limits<T>::quiet_NaN(),
      numeric_limits<T>::quiet_NaN(),
      numeric_limits<T>::quiet_NaN()
    };
  }

  auto sum = double{0};
  auto min = numeric_limits<double>::max();
  auto max = numeric_limits<double>::lowest();
  for (auto i : v) {
    sum += static_cast<double>(i);
    min = std::min(min, static_cast<double>(i));
    max = std::max(max, static_cast<double>(i));
  }
  auto avg = sum / v.size();

  auto variance = double{0};
  for (auto i : v) {
    variance += (i - avg) * (i - avg);
  }
  variance /= v.size();
  auto std_dev = sqrt(variance);

  return Stats{min, max, avg, std_dev};
}


template<typename T, std::size_t N> inline
std::vector<T> inputBuffer(
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances)
{
  using namespace std;

  if (frozen_indices.size() != frozen_distances.size()) {
    throw runtime_error("indices/distances size mismatch");
  }

  auto input_buffer = vector<T>(
    linearSize(grid_size),
    numeric_limits<T>::quiet_NaN());
  auto input_grid = Grid<T, 2>(grid_size, input_buffer.front());
  for (auto i = size_t{0}; i < frozen_indices.size(); ++i) {
    input_grid.cell(frozen_indices[i]) = frozen_distances[i];
  }
  return input_buffer;
}


template<typename T, std::size_t N> inline
std::vector<T> errorBuffer(
  std::array<std::size_t, N> const& grid_size,
  std::vector<T> const& distance_ground_truth_buffer,
  std::vector<T> const& distance_buffer)
{
  using namespace std;

  if (linearSize(grid_size) != distance_buffer.size() ||
      distance_buffer.size() != distance_ground_truth_buffer.size()) {
    throw runtime_error("distance buffers size mismatch");
  }

  auto r = vector<T>(linearSize(grid_size));
  transform(
    begin(distance_buffer), end(distance_buffer),
    begin(distance_ground_truth_buffer),
    begin(r),
    [](auto const d, auto const d_gt){ return d - d_gt; });
  return r;
}


template<typename T, std::size_t N> inline
std::vector<std::array<T, N>> distance_gradients(
  std::vector<T>& distance_buffer,
  std::array<std::size_t, N> const& grid_size,
  std::array<T, N> const& voxel_size)
{
  using namespace std;

  auto const linear_size = linearSize(grid_size);
  if (linear_size != distance_buffer.size()) {
    throw runtime_error("grid/buffer size mismatch");
  }

  Grid<T, N> distance_grid(grid_size, distance_buffer.front());

  auto grad_buffer = vector<array<T, N>>(linear_size);
  auto grad_grid = Grid<array<T, N>, N>(grid_size, grad_buffer.front());

  for (int32_t i = 0; static_cast<size_t>(i) < grid_size[0]; ++i) {
    for (int32_t j = 0; static_cast<size_t>(j) < grid_size[1]; ++j) {
      auto const index = array<int32_t, N>{{i, j}};
      grad_grid.cell(index) = gradient(
        distance_grid,
        index,
        voxel_size,
        Neighborhood<N>::offsets());
    }
  }
  return grad_buffer;
}


template<typename T> inline
void circleUnsignedFrozenCells(
  std::array<T, 2> const& center,
  T const radius,
  std::array<std::size_t, 2> const& grid_size,
  std::array<T, 2> const& voxel_size,
  std::vector<std::array<std::int32_t, 2>>* frozen_indices,
  std::vector<T>* frozen_distances,
  std::vector<T>* distance_ground_truth_buffer = nullptr)
{
  using namespace std;

  auto distance_ground_truth_grid = unique_ptr<Grid<T, 2>>();
  if (distance_ground_truth_buffer != nullptr) {
    distance_ground_truth_buffer->resize(linearSize(grid_size));
    distance_ground_truth_grid.reset(
      new Grid<T, 2>(grid_size, distance_ground_truth_buffer->front()));
  }

  for (int32_t i = 0; i < static_cast<int32_t>(grid_size[0]); ++i) {
    for (int32_t j = 0; j < static_cast<int32_t>(grid_size[1]); ++j) {
      auto const index = array<int32_t, 2>{{i, j}};
      auto const cell_corners = array<array<T, 2>, 4>
      {{
        {{(i + 0) * voxel_size[0], (j + 0) * voxel_size[1]}},
        {{(i + 1) * voxel_size[0], (j + 0) * voxel_size[1]}},
        {{(i + 0) * voxel_size[0], (j + 1) * voxel_size[1]}},
        {{(i + 1) * voxel_size[0], (j + 1) * voxel_size[1]}}
      }};

      auto inside_count = int32_t{0};
      auto outside_count = int32_t{0};

      for (auto const& cell_corner : cell_corners) {
        auto const d = distance(center, cell_corner);
        if (d < radius) {
          ++inside_count;
        }
        else {
          ++outside_count;
        }
      }

      auto const cell_center = array<T, 2>{{
        (i + T(0.5)) * voxel_size[0],
        (j + T(0.5)) * voxel_size[1]
      }};

      if (inside_count > 0 && outside_count > 0) {
        // Sphere passes through this cell. Compute and store distance
        // to cell center.
        frozen_indices->push_back(index);
        frozen_distances->push_back(fabs(distance(center, cell_center) - radius));
      }

      if (distance_ground_truth_buffer != nullptr) {
        auto const delta = sub(center, cell_center);
        distance_ground_truth_grid->cell(index) = fabs(magnitude(delta) - radius);
      }
    }
  }
}


template<typename T> inline
void circleSignedFrozenCells(
  std::array<T, 2> const& center,
  T const radius,
  std::array<std::size_t, 2> const& grid_size,
  std::array<T, 2> const& voxel_size,
  std::vector<std::array<std::int32_t, 2>>* frozen_indices,
  std::vector<T>* frozen_distances,
  std::vector<std::array<T, 2>>* normals,
  std::vector<T>* distance_ground_truth_buffer = nullptr)
{
  using namespace std;

  auto distance_ground_truth_grid = unique_ptr<Grid<T, 2>>();
  if (distance_ground_truth_buffer != nullptr) {
    distance_ground_truth_buffer->resize(linearSize(grid_size));
    distance_ground_truth_grid.reset(
      new Grid<T, 2>(grid_size, distance_ground_truth_buffer->front()));
  }

  for (int32_t i = 0; i < static_cast<int32_t>(grid_size[0]); ++i) {
    for (int32_t j = 0; j < static_cast<int32_t>(grid_size[1]); ++j) {
      auto const index = array<int32_t, 2>{{i, j}};
      auto const cell_corners = array<array<T, 2>, 4>
      {{
        {{(i + 0) * voxel_size[0], (j + 0) * voxel_size[1]}},
        {{(i + 1) * voxel_size[0], (j + 0) * voxel_size[1]}},
        {{(i + 0) * voxel_size[0], (j + 1) * voxel_size[1]}},
        {{(i + 1) * voxel_size[0], (j + 1) * voxel_size[1]}}
      }};

      auto inside_count = int32_t{0};
      auto outside_count = int32_t{0};

      for (auto const& cell_corner : cell_corners) {
        auto const d = distance(center, cell_corner);
        if (d < radius) {
          ++inside_count;
        }
        else {
          ++outside_count;
        }
      }

      auto const cell_center = array<T, 2>{{
        (i + T(0.5)) * voxel_size[0],
        (j + T(0.5)) * voxel_size[1]
      }};

      if (inside_count > 0 && outside_count > 0) {
        frozen_indices->push_back(index);
        frozen_distances->push_back(distance(cell_center, center) - radius);
        normals->push_back(normalized(sub(cell_center, center)));
      }

      if (distance_ground_truth_buffer != nullptr) {
        auto const delta = sub(center, cell_center);
        distance_ground_truth_grid->cell(index) = magnitude(delta) - radius;
      }
    }
  }
}

} // namespace detail


template<typename T, std::size_t N>
struct GradientMagnitudeStats
{
  static const std::size_t Size = N;

  double min;
  double max;
  double avg;
  double std_dev;

  std::array<std::size_t, N> grid_size;
  std::vector<T> input_buffer;
  std::vector<T> distance_buffer;
  std::vector<std::array<T, N>> grad_buffer;
};


template<typename T, std::size_t N> inline
GradientMagnitudeStats<T, N> UnsignedGradientMagnitudeStats()
{
  using namespace std;
  using namespace detail;

  auto const center = filledArray<N>(T(0.5));
  auto const radius = T(0.25);
  auto const grid_size = filledArray<N>(size_t{100});
  auto const voxel_size = filledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  circleUnsignedFrozenCells(
    center,
    radius,
    grid_size,
    voxel_size,
    &frozen_indices,
    &frozen_distances);

  auto distance_buffer = unsignedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances);

  auto const grad_buffer = distance_gradients(
    distance_buffer,
    grid_size,
    voxel_size);

  auto stats = detail::stats(transformedVector<array<T, N>, T>(
    grad_buffer, [](array<T, N> const& g) { return magnitude(g); }));

  return GradientMagnitudeStats<T, N>{
    stats.min,
    stats.max,
    stats.avg,
    stats.std_dev,
    grid_size,
    move(inputBuffer(grid_size, frozen_indices, frozen_distances)),
    move(distance_buffer),
    move(grad_buffer)};
}


template<typename T, std::size_t N> inline
GradientMagnitudeStats<T, N> SignedGradientMagnitudeStats()
{
  using namespace std;
  using namespace detail;

  auto const center = filledArray<N>(T(0.5));
  auto const radius = T(0.25);
  auto const grid_size = filledArray<N>(size_t{100});
  auto const voxel_size = filledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  circleSignedFrozenCells(
    center,
    radius,
    grid_size,
    voxel_size,
    &frozen_indices,
    &frozen_distances,
    &normals);

  auto distance_buffer = signedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances,
    normals);

  auto const grad_buffer = distance_gradients(
    distance_buffer,
    grid_size,
    voxel_size);

  auto stats = detail::stats(transformedVector<array<T, N>, T>(
    grad_buffer, [](array<T, N> const& g) { return magnitude(g); }));

  return GradientMagnitudeStats<T, N>{
    stats.min,
    stats.max,
    stats.avg,
    stats.std_dev,
    grid_size,
    move(inputBuffer(grid_size, frozen_indices, frozen_distances)),
    move(distance_buffer),
    move(grad_buffer)};
}


template<typename T, std::size_t N>
struct DistanceValueStats
{
  static const std::size_t Size = N;
  double min_error;
  double max_error;
  double avg_error;
  double std_dev_error;

  std::array<std::size_t, N> grid_size;

  std::vector<T> input_buffer;
  std::vector<T> distance_buffer;
  std::vector<T> distance_ground_truth_buffer;
  std::vector<T> error_buffer;
};


template<typename T, std::size_t N> inline
DistanceValueStats<T, N> UnsignedDistanceValueStats()
{
  using namespace std;
  using namespace detail;

  auto const center = filledArray<N>(T(0.5));
  auto const radius = T(0.25);
  auto const grid_size = filledArray<N>(size_t{100});
  auto const voxel_size = filledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto distance_ground_truth_buffer = vector<T>();
  circleUnsignedFrozenCells(
    center,
    radius,
    grid_size,
    voxel_size,
    &frozen_indices,
    &frozen_distances,
    &distance_ground_truth_buffer);

  auto const input_buffer = inputBuffer(
    grid_size,
    frozen_indices,
    frozen_distances);

  auto distance_buffer = unsignedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances);

  auto const error_buffer = errorBuffer(
    grid_size,
    distance_buffer,
    distance_ground_truth_buffer);

  auto const stats = detail::stats(error_buffer);

  return DistanceValueStats<T, N>{
    stats.min,
    stats.max,
    stats.avg,
    stats.std_dev,
    grid_size,
    input_buffer,
    distance_buffer,
    distance_ground_truth_buffer,
    error_buffer};
}


template<typename T, std::size_t N> inline
DistanceValueStats<T, N> SignedDistanceValueStats()
{
  using namespace std;
  using namespace detail;

  auto const center = filledArray<N>(T(0.5));
  auto const radius = T(0.25);
  auto const grid_size = filledArray<N>(size_t{100});
  auto const voxel_size = filledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  auto distance_ground_truth_buffer = vector<T>();
  circleSignedFrozenCells(
    center,
    radius,
    grid_size,
    voxel_size,
    &frozen_indices,
    &frozen_distances,
    &normals,
    &distance_ground_truth_buffer);

  auto distance_buffer = signedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances,
    normals);

  auto const input_buffer = inputBuffer(
    grid_size,
    frozen_indices,
    frozen_distances);

  auto const error_buffer = errorBuffer(
    grid_size,
    distance_buffer,
    distance_ground_truth_buffer);

  auto const stats = detail::stats(error_buffer);

  return DistanceValueStats<T, N>{
    stats.min,
    stats.max,
    stats.avg,
    stats.std_dev,
    grid_size,
    input_buffer,
    distance_buffer,
    distance_ground_truth_buffer,
    error_buffer};
}

} // namespace test
} // namespace fmm
} // namespace thinks

#endif // THINKS_TESTFASTMARCHINGMETHOD_HPP_INCLUDED
