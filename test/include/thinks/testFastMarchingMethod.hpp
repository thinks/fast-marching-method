#ifndef THINKS_TESTFASTMARCHINGMETHOD_HPP_INCLUDED
#define THINKS_TESTFASTMARCHINGMETHOD_HPP_INCLUDED

#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <functional>
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
  cell_type* const cells_;
};

namespace detail {

//! Returns base^exponent as a compile-time constant.
constexpr std::size_t static_pow(std::size_t base, std::size_t const exponent)
{
  using namespace std;

  // NOTE: Cannot use loops in constexpr functions in C++11, have to use
  // recursion here.
  return exponent == size_t{0} ? size_t{1} : base * static_pow(base, exponent - 1);
}


//! Returns the product of the elements in array @a a.
//! Note: Not checking for overflow here!
template<std::size_t N> inline
std::size_t LinearSize(std::array<std::size_t, N> const& a)
{
  using namespace std;
  return accumulate(begin(a), end(a), size_t{1}, multiplies<size_t>());
}


template<typename T, typename R, typename U> inline
std::vector<R> TransformedVector(std::vector<T> const& v, U const unary_op)
{
  using namespace std;

  auto r = vector<R>(v.size());
  transform(begin(v), end(v), begin(r), unary_op);
  return r;
}


template<std::size_t N, typename T> inline
std::array<T, N> FilledArray(T const v)
{
  using namespace std;

  auto r = array<T, N>{};
  fill(begin(r), end(r), v);
  return r;
}


template<typename T, std::size_t N> inline
std::array<T, N> add(std::array<T, N> const& u, std::array<T, N> const& v)
{
  using namespace std;

  auto r = array<T, N>{};
  transform(begin(u), end(u), begin(v), begin(r),
            [](auto const ux, auto const vx){ return ux + vx; });
  return r;
}


template<typename T, std::size_t N> inline
std::array<T, N> sub(std::array<T, N> const& u, std::array<T, N> const& v)
{
  using namespace std;

  auto r = array<T, N>{};
  transform(begin(u), end(u), begin(v), begin(r),
            [](auto const ux, auto const vx){ return ux - vx; });
  return r;
}


template<typename T, std::size_t N> inline
std::array<T, N> mult(T const t, std::array<T, N> const& v)
{
  using namespace std;

  auto r = v;
  for (auto i = size_t{0}; i < N; ++i) {
    v *= t;
  }
  return r;
}


template<typename T, std::size_t N> inline
T dot(std::array<T, N> const& u, std::array<T, N> const& v)
{
  using namespace std;

  auto sum = T(0);
  for (auto i = size_t{0}; i < N; ++i) {
    sum += u[i] * v[i];
  }
  return sum;
}


template<typename T> inline
std::array<T, 2> solveQuadratic(T const a, T const b, T const c)
{
  using namespace std;

  T const eps = 1e-9;

  assert(abs(a) > eps);

  T const discr = b * b - T(4) * a * c;
  if (discr < T(0)) {
    return {{numeric_limits<T>::quiet_NaN(), numeric_limits<T>::quiet_NaN()}};
  }
  else if (discr < eps) {
    T const r = T(-0.5) * b / a;
    return {{r, r}};
  }

  T const q = (b > 0) ?
    T(-0.5) * (b + sqrt(discr)) : T(-0.5) * (b - sqrt(discr));
  T const x0 = q / a;
  T const x1 = c / q;
  return {{min(x0, x1), max(x0, x1)}};
}


template<typename T, std::size_t N> inline
std::array<T, N> rayHyperSphereIntersection(
  std::array<T, N> const& rayOrigin,
  std::array<T, N> const& rayDirection,
  std::array<T, N> const& center,
  T const radius)
{
  using namespace std;

  auto const p = sub(rayOrigin, center);
  auto const a = dot(rayDirection, rayDirection);
  auto const b = T(2) * dot(rayDirection, p);
  auto const c = dot(p, p) - radius * radius;

  auto const q = solveQuadratic(a, b, c);

  if (any_of(begin(q), end(q), [](auto const x) { return isnan(x); })) {
    // No intersection.
    return FilledArray<T, N>(numeric_limits<T>::quiet_NaN());
  }

  if (q[0] < T(0)) {
    if (q[1] < T(0)) {
      // Both intersections behind ray origin.
      return FilledArray<T, N>(numeric_limits<T>::quiet_NaN());
    }

    return add(rayOrigin, mult(q[1], rayDirection));
  }

  return add(rayOrigin, mult(q[0], rayDirection));
}


template<std::size_t N>
struct Neighborhood
{
  static std::array<std::array<std::int32_t, N>, 2 * N> offsets()
  {
    using namespace std;

    array<array<int32_t, N>, 2 * N> n;
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


template<std::size_t N>
class IndexIterator
{
public:
  explicit IndexIterator(std::array<std::size_t, N> const& size)
    : size_(size)
  {
    using namespace std;

    for (auto const s : size) {
      if (s < 1) {
        throw runtime_error("zero size element");
      }
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

    // Use central difference if possible.
    auto const pos_inside =
      0 <= pos_index[i] && static_cast<size_t>(pos_index[i]) < size[i];
    auto const neg_inside =
      0 <= neg_index[i] && static_cast<size_t>(neg_index[i]) < size[i];
    if (pos_inside && neg_inside) {
      grad[i] = (grid.cell(pos_index) - grid.cell(neg_index)) / (T(2) * dx[i]);
    }
    else if (pos_inside) {
      grad[i] = (grid.cell(pos_index) - grid.cell(index)) / dx[i];
    }
    else if (neg_inside) {
      grad[i] = (grid.cell(index) - grid.cell(neg_index)) / dx[i];
    }
    else {
      grad[i] = numeric_limits<T>::quiet_NaN();
    }
  }
  return grad;
}


template<typename T, std::size_t N> inline
T Magnitude(std::array<T, N> const& v)
{
  using namespace std;

  auto mag_squared = T{0};
  for (auto i = size_t{0}; i < N; ++i) {
    mag_squared += v[i] * v[i];
  }
  return sqrt(mag_squared);
}


template<typename T, std::size_t N> inline
T Distance(std::array<T, N> const& u, std::array<T, N> const& v)
{
  using namespace std;

  auto distance_squared = T{0};
  for (auto i = size_t{0}; i < N; ++i) {
    auto const delta = u[i] - v[i];
    distance_squared += delta * delta;
  }
  return sqrt(distance_squared);
}


template<typename T, std::size_t N> inline
std::array<T, N> Normalized(std::array<T, N> const& v)
{
  using namespace std;

  auto n = v;
  auto const mag_factor = T{1} / Magnitude(v);
  for (auto i = size_t{0}; i < N; ++i) {
    n[i] *= mag_factor;
  }
  return n;
}


struct ErrorStats
{
  double min_abs_error;
  double max_abs_error;
  double avg_abs_error;
  double std_dev_abs_error;
};


template<typename T> inline
ErrorStats ErrorStatistics(std::vector<T> const& errors)
{
  using namespace std;

  static_assert(is_floating_point<T>::value, "value type must be floating point");

  if (errors.empty()) {
    return ErrorStats{
      numeric_limits<T>::quiet_NaN(),
      numeric_limits<T>::quiet_NaN(),
      numeric_limits<T>::quiet_NaN(),
      numeric_limits<T>::quiet_NaN()
    };
  }

  auto sum_abs_error = double{0};
  auto min_abs_error = numeric_limits<double>::max();
  auto max_abs_error = double{0};
  for (auto e : errors) {
    auto const ae = static_cast<double>(fabs(e));
    sum_abs_error += ae;
    min_abs_error = min(min_abs_error, ae);
    max_abs_error = max(max_abs_error, ae);
  }
  auto avg = sum_abs_error / errors.size();

  auto variance = double{0};
  for (auto e : errors) {
    variance += (fabs(e) - avg) * (fabs(e) - avg);
  }
  variance /= errors.size();
  auto std_dev = sqrt(variance);

  return ErrorStats{min_abs_error, max_abs_error, avg, std_dev};
}


template<typename T, std::size_t N> inline
std::vector<T> InputBuffer(
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>> const& frozen_indices,
  std::vector<T> const& frozen_distances)
{
  using namespace std;

  if (frozen_indices.size() != frozen_distances.size()) {
    throw runtime_error("indices/distances size mismatch");
  }

  auto input_buffer = vector<T>(LinearSize(grid_size), numeric_limits<T>::max());
  auto input_grid = Grid<T, N>(grid_size, input_buffer.front());
  for (auto i = size_t{0}; i < frozen_indices.size(); ++i) {
    input_grid.cell(frozen_indices[i]) = frozen_distances[i];
  }
  return input_buffer;
}


template<typename T> inline
std::vector<T> ErrorBuffer(
  std::vector<T> const& ground_truth_buffer,
  std::vector<T> const& values_buffer)
{
  using namespace std;

  if (ground_truth_buffer.size() != values_buffer.size()) {
    throw runtime_error("buffer size mismatch");
  }

  auto r = vector<T>(values_buffer.size());
  transform(
    begin(values_buffer), end(values_buffer),
    begin(ground_truth_buffer),
    begin(r),
    [](auto const v, auto const gt){ return v - gt; });
  return r;
}


template<typename T, std::size_t N> inline
std::vector<std::array<T, N>> DistanceGradients(
  std::vector<T>& distance_buffer,
  std::array<std::size_t, N> const& grid_size,
  std::array<T, N> const& voxel_size)
{
  using namespace std;

  auto const linear_size = LinearSize(grid_size);
  if (linear_size != distance_buffer.size()) {
    throw runtime_error("grid/buffer size mismatch");
  }

  Grid<T, N> distance_grid(grid_size, distance_buffer.front());

  auto grad_buffer = vector<array<T, N>>(linear_size);
  auto grad_grid = Grid<array<T, N>, N>(grid_size, grad_buffer.front());

  auto index_iter = IndexIterator<N>(grid_size);
  auto valid_index = bool{true};
  while (valid_index) {
    auto const index = index_iter.index();
    grad_grid.cell(index) = gradient(
      distance_grid,
      index,
      voxel_size,
      Neighborhood<N>::offsets());
    valid_index = index_iter.Next();
  }

  return grad_buffer;
}


template<typename T, std::size_t N> inline
std::array<T, N> CellCenter(
  std::array<std::int32_t, N> const& index,
  std::array<T, N> const& voxel_size)
{
  using namespace std;

  auto cell_center = array<T, N>{};
  for (auto i = size_t{0}; i < N; ++i) {
    cell_center[i] = (index[i] + T(0.5)) * voxel_size[i];
  }
  return cell_center;
}


template<typename T, std::size_t N> inline
std::array<std::array<T, N>, static_pow(2, N)> CellCorners(
  std::array<std::int32_t, N> const& index,
  std::array<T, N> const& voxel_size)
{
  using namespace std;

  auto cell_corners = array<array<T, N>, static_pow(2, N)>{};

  for (auto i = size_t{0}; i < static_pow(2, N); ++i) {
    auto const bits = bitset<N>(i);
    for (auto k = size_t{0}; k < N; ++k) {
      cell_corners[i][k] =
        (index[k] + static_cast<int32_t>(bits[k])) * voxel_size[k];
    }
  }
  return cell_corners;
}


template<typename T, std::size_t N> inline
void HyperSphereFrozenCells(
  std::array<T, N> const& center,
  T const radius,
  std::array<std::size_t, N> const& grid_size,
  std::array<T, N> const& voxel_size,
  std::vector<std::array<std::int32_t, N>>* frozen_indices,
  std::vector<T>* frozen_distances,
  std::vector<std::array<T, N>>* normals,
  std::vector<T>* distance_ground_truth_buffer = nullptr,
  std::function<T(T)> const ground_truth_distance_op =
    [](T const d) { return d; })
{
  using namespace std;

  auto distance_ground_truth_grid = unique_ptr<Grid<T, N>>();
  if (distance_ground_truth_buffer != nullptr) {
    distance_ground_truth_buffer->resize(LinearSize(grid_size));
    distance_ground_truth_grid.reset(
      new Grid<T, N>(grid_size, distance_ground_truth_buffer->front()));
  }

  auto index_iter = IndexIterator<N>(grid_size);
  auto valid_index = bool{true};
  while (valid_index) {
    auto const index = index_iter.index();
    auto const cell_corners = CellCorners(index, voxel_size);

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

    auto const cell_center = CellCenter(index, voxel_size);
    auto const cell_distance = Distance(center, cell_center) - radius;

    if (inside_count > 0 && outside_count > 0) {
      // The inferface passes through this cell.
      frozen_indices->push_back(index);
      frozen_distances->push_back(cell_distance);
      normals->push_back(Normalized(sub(cell_center, center)));
    }

    if (distance_ground_truth_buffer != nullptr) {
      distance_ground_truth_grid->cell(index) =
        ground_truth_distance_op(cell_distance);
    }

    valid_index = index_iter.Next();
  }
}

} // namespace detail


template<typename T, std::size_t N>
struct GradientMagnitudeStats
{
  static const std::size_t kSize = N;

  double min_abs_error;
  double max_abs_error;
  double avg_abs_error;
  double std_dev_abs_error;

  double duration_in_s;

  std::array<std::size_t, N> grid_size;
  std::vector<T> input_buffer;
  std::vector<T> distance_buffer;
  std::vector<std::array<T, N>> grad_buffer;
  std::vector<T> error_buffer;
};


template<typename T, std::size_t N> inline
GradientMagnitudeStats<T, N> UnsignedGradientMagnitudeStats()
{
  using namespace std;
  using namespace detail;

  auto const center = FilledArray<N>(T(0.5));
  auto const radius = T(0.25);
  auto const grid_size = FilledArray<N>(size_t{100});
  auto const voxel_size = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  HyperSphereFrozenCells<T, N>(
    center,
    radius,
    grid_size,
    voxel_size,
    &frozen_indices,
    &frozen_distances,
    &normals);

  // TMP!!
  for (auto iter = begin(frozen_distances); iter != end(frozen_distances); ++iter) {
    *iter = fabs(*iter);
  }

  auto const start_time = chrono::system_clock::now();
  auto distance_buffer = UnsignedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances);
  auto const end_time = chrono::system_clock::now();

  auto const input_buffer = InputBuffer(
    grid_size,
    frozen_indices,
    frozen_distances);

  auto const grad_buffer = DistanceGradients(
    distance_buffer,
    grid_size,
    voxel_size);

  auto ground_truth_buffer = vector<T>(grad_buffer.size());
  fill(begin(ground_truth_buffer), end(ground_truth_buffer), T(1));
  auto const error_buffer = ErrorBuffer(
    TransformedVector<array<T, N>, T>(
      grad_buffer,
      [](array<T, N> const& g) { return Magnitude(g); }),
    ground_truth_buffer);
  auto stats = ErrorStatistics(error_buffer);

  return GradientMagnitudeStats<T, N>{
    stats.min_abs_error,
    stats.max_abs_error,
    stats.avg_abs_error,
    stats.std_dev_abs_error,
    chrono::duration<double>(end_time - start_time).count(),
    grid_size,
    input_buffer,
    distance_buffer,
    grad_buffer,
    error_buffer};
}


template<typename T, std::size_t N> inline
GradientMagnitudeStats<T, N> SignedGradientMagnitudeStats()
{
  using namespace std;
  using namespace detail;

  auto const center = FilledArray<N>(T(0.5));
  auto const radius = T(0.25);
  auto const grid_size = FilledArray<N>(size_t{100});
  auto const voxel_size = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  HyperSphereFrozenCells<T, N>(
    center,
    radius,
    grid_size,
    voxel_size,
    &frozen_indices,
    &frozen_distances,
    &normals);

  auto const start_time = chrono::system_clock::now();
  auto distance_buffer = SignedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances);
  auto const end_time = chrono::system_clock::now();

  auto const input_buffer = InputBuffer(
    grid_size,
    frozen_indices,
    frozen_distances);

  auto const grad_buffer = DistanceGradients(
    distance_buffer,
    grid_size,
    voxel_size);

  auto ground_truth_buffer = vector<T>(grad_buffer.size());
  fill(begin(ground_truth_buffer), end(ground_truth_buffer), T(1));
  auto const error_buffer = ErrorBuffer(
    TransformedVector<array<T, N>, T>(
      grad_buffer,
      [](array<T, N> const& g) { return Magnitude(g); }),
    ground_truth_buffer);
  auto stats = ErrorStatistics(error_buffer);

  return GradientMagnitudeStats<T, N>{
    stats.min_abs_error,
    stats.max_abs_error,
    stats.avg_abs_error,
    stats.std_dev_abs_error,
    chrono::duration<double>(end_time - start_time).count(),
    grid_size,
    input_buffer,
    distance_buffer,
    grad_buffer,
    error_buffer};
}


template<typename T, std::size_t N>
struct DistanceValueStats
{
  static const std::size_t kSize = N;

  double min_abs_error;
  double max_abs_error;
  double avg_abs_error;
  double std_dev_abs_error;

  double duration_in_s;

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

  auto const center = FilledArray<N>(T(0.5));
  auto const radius = T(0.25);
  auto const grid_size = FilledArray<N>(size_t{100});
  auto const voxel_size = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  auto distance_ground_truth_buffer = vector<T>();
  HyperSphereFrozenCells<T, N>(
    center,
    radius,
    grid_size,
    voxel_size,
    &frozen_indices,
    &frozen_distances,
    &normals,
    &distance_ground_truth_buffer,
    [](auto const d) { return fabs(d); });

  // TMP!!
  for (auto iter = begin(frozen_distances); iter != end(frozen_distances); ++iter) {
    *iter = fabs(*iter);
  }

  auto const start_time = chrono::system_clock::now();
  auto distance_buffer = UnsignedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances);
  auto const end_time = chrono::system_clock::now();

  auto const input_buffer = InputBuffer(
    grid_size,
    frozen_indices,
    frozen_distances);

  auto const error_buffer = ErrorBuffer(
    distance_buffer,
    distance_ground_truth_buffer);

  auto const stats = ErrorStatistics(error_buffer);

  return DistanceValueStats<T, N>{
    stats.min_abs_error,
    stats.max_abs_error,
    stats.avg_abs_error,
    stats.std_dev_abs_error,
    chrono::duration<double>(end_time - start_time).count(),
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

  auto const center = FilledArray<N>(T(0.5));
  auto const radius = T(0.25);
  auto const grid_size = FilledArray<N>(size_t{100});
  auto const voxel_size = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  auto distance_ground_truth_buffer = vector<T>();
  HyperSphereFrozenCells<T, N>(
    center,
    radius,
    grid_size,
    voxel_size,
    &frozen_indices,
    &frozen_distances,
    &normals,
    &distance_ground_truth_buffer);

  auto const start_time = chrono::system_clock::now();
  auto distance_buffer = SignedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances);
  auto const end_time = chrono::system_clock::now();

  auto const input_buffer = InputBuffer(
    grid_size,
    frozen_indices,
    frozen_distances);

  auto const error_buffer = ErrorBuffer(
    distance_buffer,
    distance_ground_truth_buffer);

  auto const stats = ErrorStatistics(error_buffer);

  return DistanceValueStats<T, N>{
    stats.min_abs_error,
    stats.max_abs_error,
    stats.avg_abs_error,
    stats.std_dev_abs_error,
    chrono::duration<double>(end_time - start_time).count(),
    grid_size,
    input_buffer,
    distance_buffer,
    distance_ground_truth_buffer,
    error_buffer};
}


#if 0
template<typename T, std::size_t N> inline
DistanceValueStats<T, N> unsignedNegativeCenter()
{
  using namespace std;
  using namespace detail;

  auto const grid_size = FilledArray<N>(size_t{100});
  auto const voxel_size = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  frozen_indices.push_back({{0, 1}});
  frozen_indices.push_back({{2, 1}});
  frozen_indices.push_back({{1, 0}});
  frozen_indices.push_back({{1, 2}});
  auto frozen_distances = vector<T>();
  frozen_distances.push_back(T(0));
  frozen_distances.push_back(T(100));
  frozen_distances.push_back(T(0));
  frozen_distances.push_back(T(100));
  auto distance_buffer = UnsignedDistance<T, N>(
    grid_size,
    voxel_size,
    speed,
    frozen_indices,
    frozen_distances);

  return DistanceValueStats<T, N>{
    numeric_limits<double>::quiet_NaN(),
    numeric_limits<double>::quiet_NaN(),
    numeric_limits<double>::quiet_NaN(),
    numeric_limits<double>::quiet_NaN(),
    numeric_limits<double>::quiet_NaN(),
    grid_size,
    distance_buffer, // HACK
    distance_buffer,
    distance_buffer, // HACK
    distance_buffer};// HACK
}
#endif

} // namespace test
} // namespace fmm
} // namespace thinks

#endif // THINKS_TESTFASTMARCHINGMETHOD_HPP_INCLUDED
