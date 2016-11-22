#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <thinks/fastMarchingMethod.hpp>
#include <../../ppm-io/include/thinks/ppm.hpp> // TMP!!!
#include <iostream>


namespace {
namespace util {

#if 1 // TMP!!!
template<typename R, typename T> inline
R Clamp(T const low, T const high, T const value)
{
  using namespace std;

  return static_cast<R>(min<T>(high, max<T>(low, value)));
}

template<typename InIter, typename C> inline
std::vector<std::uint8_t> PixelsFromValues(
  InIter const in_begin,
  InIter const in_end,
  C const pixel_from_value)
{
  using namespace std;

  auto pixels = vector<uint8_t>();
  for (auto iter = in_begin; iter != in_end; ++iter) {
    auto const value = *iter;
    auto const pixel = pixel_from_value(value);
    pixels.insert(end(pixels), begin(pixel), end(pixel));
  }
  return pixels;
}

template<typename T, typename InIter> inline
std::vector<T> SignedNormalized(InIter const in_begin, InIter const in_end)
{
  using namespace std;

  auto max_pos_value = numeric_limits<T>::lowest();
  auto min_neg_value = numeric_limits<T>::max();
  for (auto iter = in_begin; iter != in_end; ++iter) {
    auto const value = *iter;
    if (value == numeric_limits<T>::max()) {
      continue;
    }

    if (value > T(0)) {
      max_pos_value = max(max_pos_value, value);
    }

    if (value < T(0)) {
      min_neg_value = min(min_neg_value, value);
    }
  }

  auto const pos_factor = max_pos_value > T(0) ? T(1) / max_pos_value : T(-1);
  auto const neg_factor = min_neg_value < T(0) ? T(-1) / min_neg_value : T(-1);

  auto r = vector<T>{};
  transform(
    in_begin,
    in_end,
    back_inserter(r),
    [=](auto const value) {
      if (value == numeric_limits<T>::max()) {
        return value;
      }

      if (value > T(0) && pos_factor > T(0)) {
        return pos_factor * value;
      }

      if (value < T(0) && neg_factor > T(0)) {
        return neg_factor * value;
      }

      return T(0);
    });
  return r;
}

template<typename T>
void writeRgbImage(
  std::string const& filename,
  std::size_t width,
  std::size_t height,
  std::vector<T> const& pixel_data)
{
  using namespace std;

  // Negative values in shades of blue, positive values in shades of red.
  // Very large values as grey.
  auto const pixel_from_value = [](T const x) {
    if (x == numeric_limits<T>::max()) {
      return array<uint8_t, 3>{{128, 128, 128}};
    }
    return x < T{0} ?
      array<uint8_t, 3>{{
        uint8_t{0},
        uint8_t{0},
        Clamp<uint8_t>(
          T{0},
          T(numeric_limits<uint8_t>::max()),
          numeric_limits<uint8_t>::max() * fabs(x))}} :
      array<uint8_t, 3>{{
        Clamp<uint8_t>(
          T{0},
          T(numeric_limits<uint8_t>::max()),
          numeric_limits<uint8_t>::max() * x),
        uint8_t{0},
        uint8_t{0}}};
  };

  auto const norm_pixel_data =
    SignedNormalized<T>(begin(pixel_data), end(pixel_data));
  thinks::ppm::writeRgbImage(
    filename,
    width,
    height,
    PixelsFromValues(
      begin(norm_pixel_data),
      end(norm_pixel_data),
      pixel_from_value));
}

#endif // TMP!!!


//! Returns base^exponent as a compile-time constant.
constexpr std::size_t static_pow(std::size_t base, std::size_t const exponent)
{
  using namespace std;

  // NOTE: Cannot use loops in constexpr functions in C++11, have to use
  // recursion here.
  return exponent == size_t{0} ? size_t{1} : base * static_pow(base, exponent - 1);
}


//! Returns the product of the elements in array @a a.
//! Note: Not checking for integer overflow here!
template<std::size_t N> inline
std::size_t LinearSize(std::array<std::size_t, N> const& a)
{
  using namespace std;
  return accumulate(begin(a), end(a), size_t{1}, multiplies<size_t>());
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


//! Returns an array with all elements initialized to have value @a v.
template<std::size_t N, typename T> inline
std::array<T, N> FilledArray(T const v)
{
  using namespace std;

  auto r = array<T, N>{};
  fill(begin(r), end(r), v);
  return r;
}


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


//! Returns a pair where the first element is true if the provided function
//! (@a func) threw an exception. The second element is set to the exception
//! message.
template<typename E, typename F>
std::pair<bool, std::string> FunctionThrows(F const func)
{
  using namespace std;

  auto thrown = false;
  auto reason = string();
  try {
    func();
  }
  catch (exception& ex)
  {
    thrown = true;
    auto const typed_ex = dynamic_cast<E*>(&ex);
    if (typed_ex == nullptr) {
      auto ss = stringstream();
      ss << "incorrect exception type: '" << typeid(ex).name() << "' "
         << "(expected '" << typeid(E).name() << "')";
      reason = ss.str();
    }
    else {
      reason = ex.what();
    }
  }

  return {thrown, reason};
}


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
//!   [0, 1]
//!   [1, 0]
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

    // Nothing was incremented.
    has_next_ = false;
    return false;
  }

private:
  std::array<std::size_t, N> const size_;
  std::array<std::int32_t, N> index_;
  bool has_next_;
};


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
template<typename T, std::size_t N, typename D> inline
void HyperSphereBoundaryCells(
  std::array<T, N> const& center,
  T const radius,
  std::array<std::size_t, N> const& grid_size,
  std::array<T, N> const& grid_spacing,
  D const distance_modifier,
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
    }

    if (distance_ground_truth_grid != nullptr) {
      distance_ground_truth_grid->Cell(index) = cell_distance;
    }

    index_iter.Next();
  }
}


//! Returns a grid buffer where boundary distances have been set. Cells that
//! are not on the boundary have the value std::numeric_limits<T>::max().
template<typename T, std::size_t N>
std::vector<T> InputBuffer(
  std::array<std::size_t, N> const& grid_size,
  std::vector<std::array<std::int32_t, N>> const& boundary_indices,
  std::vector<T> const& boundary_distances)
{
  using namespace std;

  if (boundary_indices.size() != boundary_distances.size()) {
    throw runtime_error("boundary indices/distances size mismatch");
  }

  auto input_buffer =
      vector<T>(LinearSize(grid_size), numeric_limits<T>::max());
  auto input_grid = Grid<T, N>(grid_size, input_buffer.front());
  for (auto i = size_t{0}; i < boundary_indices.size(); ++i) {
    auto const boundary_index = boundary_indices[i];
    if (!util::Inside(boundary_index, grid_size)) {
      throw runtime_error("index outside grid");
    }

    input_grid.Cell(boundary_index) = boundary_distances[i];
  }

  return input_buffer;
}


//! Returns a vector of element-wise errors. If a value is larger than the
//! corresponding ground truth a positive error is reported, if a value is
//! smaller than its corresponding ground truth a negative value is reported.
//! The unit of the errors is voxel units. For non-uniform grid spacing the
//! smallest component is used.
template<typename T, std::size_t N>
std::vector<T> ErrorBuffer(
  std::vector<T> const& value_buffer,
  std::vector<T> const& ground_truth_buffer,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  if (ground_truth_buffer.size() != value_buffer.size()) {
    throw runtime_error("buffer size mismatch");
  }

  auto const min_grid_spacing = *min_element(begin(grid_spacing), end(grid_spacing));

  auto r = vector<T>(value_buffer.size());
  transform(
    begin(value_buffer), end(value_buffer),
    begin(ground_truth_buffer),
    begin(r),
    [=](auto const v, auto const gt){ return (v - gt) / min_grid_spacing; });
  return r;
}


struct ErrorStats
{
  double min_abs_error;
  double max_abs_error;
  double avg_abs_error;
  double std_dev_abs_error;
};


//! Returns basic error metrics from the (non-empty) list of scalar error
//! values.
template<typename T> inline
ErrorStats ErrorStatistics(std::vector<T> const& errors)
{
  using namespace std;

  static_assert(is_floating_point<T>::value, "value type must be floating point");

  if (errors.empty()) {
    throw runtime_error("empty error vector");
  }

  auto sum_abs_error = double{0};
  auto min_abs_error = numeric_limits<double>::max();
  auto max_abs_error = double{0};
  for (auto const error : errors) {
    auto const abs_error = static_cast<double>(fabs(error));
    sum_abs_error += abs_error;
    min_abs_error = min(min_abs_error, abs_error);
    max_abs_error = max(max_abs_error, abs_error);
  }
  auto const avg_abs_error = sum_abs_error / errors.size();

  auto error_variance = double{0};
  for (auto const error : errors) {
    auto const delta_error = fabs(error) - avg_abs_error;
    error_variance += delta_error * delta_error;
  }
  error_variance /= errors.size();
  auto error_std_dev = sqrt(error_variance);

  return ErrorStats{min_abs_error, max_abs_error, avg_abs_error, error_std_dev};
}


template<typename T, std::size_t N>
struct DistanceStats
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


template<typename T, std::size_t N,
         typename DistanceFunction, typename DistanceModifier>
DistanceStats<T, N> HyperSphereDistanceStats(
  std::array<T, N> const sphere_center,
  T const sphere_radius,
  std::array<std::size_t, N> const grid_size,
  std::array<T, N> const grid_spacing,
  DistanceModifier const distance_modifier,
  DistanceFunction const distance_function)
{
  using namespace std;

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto distance_ground_truth_buffer = vector<T>();

  HyperSphereBoundaryCells<T, N>(
    sphere_center,
    sphere_radius,
    grid_size,
    grid_spacing,
    distance_modifier,
    &frozen_indices,
    &frozen_distances,
    &distance_ground_truth_buffer);

  auto const start_time = chrono::system_clock::now();
  auto distance_buffer =
    distance_function(
      frozen_indices,
      frozen_distances);
  auto const end_time = chrono::system_clock::now();

  auto const input_buffer = InputBuffer(
    grid_size,
    frozen_indices,
    frozen_distances);

  auto const error_buffer = ErrorBuffer(
    distance_buffer,
    distance_ground_truth_buffer,
    grid_spacing);

  auto const error_stats = ErrorStatistics(error_buffer);

  return DistanceStats<T, N>{
    error_stats.min_abs_error,
    error_stats.max_abs_error,
    error_stats.avg_abs_error,
    error_stats.std_dev_abs_error,
    chrono::duration<double>(end_time - start_time).count(),
    grid_size,
    input_buffer,
    distance_buffer,
    distance_ground_truth_buffer,
    error_buffer};
}


template<typename T, std::size_t N>
std::array<T, N> Gradient(
  std::array<std::int32_t, N> const& index,
  Grid<T, N> const& grid,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  static_assert(is_floating_point<T>::value,
                "scalar type must be floating point");

  array<T, N> grad;
  auto const size = grid.size();

  for (size_t i = 0; i < N; ++i) {
    auto pos_index = index;
    auto neg_index = index;
    pos_index[i] += size_t{1};
    neg_index[i] -= size_t{1};

    auto const cell = grid.Cell(index);

    // Use central difference if possible.
    auto const pos_inside =
      0 <= pos_index[i] && static_cast<size_t>(pos_index[i]) < size[i];
    auto const neg_inside =
      0 <= neg_index[i] && static_cast<size_t>(neg_index[i]) < size[i];
    if (pos_inside && neg_inside) {
      auto const pos_cell = grid.Cell(pos_index);
      auto const neg_cell = grid.Cell(neg_index);
      if (!(pos_cell > cell && neg_cell > cell)) {
        grad[i] = (pos_cell - neg_cell) / (T{2} * grid_spacing[i]);
      }
      else {
        grad[i] = (pos_cell - cell) / grid_spacing[i];
      }
    }
    else if (pos_inside) {
      auto const pos_cell = grid.Cell(pos_index);
      grad[i] = (pos_cell - cell) / grid_spacing[i];
    }
    else if (neg_inside) {
      auto const neg_cell = grid.Cell(neg_index);
      grad[i] = (cell - neg_cell) / grid_spacing[i];
    }
    else {
      throw runtime_error("invalid gradient");
    }
  }
  return grad;
}


template<typename T, std::size_t N> inline
std::vector<std::array<T, N>> DistanceGradients(
  Grid<T, N> const& distance_grid,
  std::array<T, N> const& grid_spacing)
{
  using namespace std;

  auto grad_buffer = vector<array<T, N>>(LinearSize(distance_grid.size()));
  auto grad_grid = Grid<array<T, N>, N>(distance_grid.size(), grad_buffer.front());

  auto index_iter = IndexIterator<N>(grad_grid.size());
  while (index_iter.has_next()) {
    auto const index = index_iter.index();
    grad_grid.Cell(index) = Gradient(index, distance_grid, grid_spacing);
    index_iter.Next();
  }

  return grad_buffer;
}


template<typename T, std::size_t N>
std::vector<T> GradientMagnitudeErrors(
  Grid<T, N> const& grad_mag_grid,
  std::vector<std::array<int32_t, N>> const& boundary_indices)
{
  using namespace std;

  // Create a mask where 1 means that the cell is on the boundary and
  // 0 means that the cell is not on the boundary.
  auto boundary_buffer =
    vector<uint8_t>(LinearSize(grad_mag_grid.size()), uint8_t{0});
  auto boundary_grid =
    Grid<uint8_t, N>(grad_mag_grid.size(), boundary_buffer.front());
  for (auto const& boundary_index : boundary_indices) {
    boundary_grid.Cell(boundary_index) = uint8_t{1};
  }

  auto grad_mag_error_buffer = vector<T>(LinearSize(grad_mag_grid.size()), T{0});
  auto grad_mag_error_grid =
    Grid<T, N>(grad_mag_grid.size(), grad_mag_error_buffer.front());

  auto index_iter = IndexIterator<N>(grad_mag_error_grid.size());
  while (index_iter.has_next()) {
    auto const index = index_iter.index();
    // Don't compute errors for boundary cells.
    if (boundary_grid.Cell(index) == uint8_t{0}) {
      grad_mag_error_grid.Cell(index) = grad_mag_grid.Cell(index) - T{1};
    }
    index_iter.Next();
  }

  return grad_mag_error_buffer;
}


template<typename T, std::size_t N>
struct UnsignedDistanceAccuracyThresholds;

template<>
struct UnsignedDistanceAccuracyThresholds<float, 2>
{
  static constexpr double min_abs_error_threshold() { return double{0}; }
  static constexpr double max_abs_error_threshold() { return double{1}; }
  static constexpr double avg_abs_error_threshold() { return double{1}; }
  static constexpr double std_dev_abs_error_threshold() { return double{1}; }
};

template<>
struct UnsignedDistanceAccuracyThresholds<double, 2>
{
  static constexpr double min_abs_error_threshold() { return double{0}; }
  static constexpr double max_abs_error_threshold() { return double{1}; }
  static constexpr double avg_abs_error_threshold() { return double{1}; }
  static constexpr double std_dev_abs_error_threshold() { return double{1}; }
};

} // namespace util

template<typename S, std::size_t N>
struct ScalarDimensionPair
{
  typedef S ScalarType;
  static constexpr std::size_t kDimension = N;
};

typedef ::testing::Types<
  ScalarDimensionPair<float, 1>,
  ScalarDimensionPair<float, 2>,
  ScalarDimensionPair<float, 3>,
  ScalarDimensionPair<float, 4>,
  ScalarDimensionPair<double, 1>,
  ScalarDimensionPair<double, 2>,
  ScalarDimensionPair<double, 3>,
  ScalarDimensionPair<double, 4>> EikonalSolverTypes;

typedef ::testing::Types<
  ScalarDimensionPair<float, 1>,
  ScalarDimensionPair<float, 2>,
  ScalarDimensionPair<float, 3>,
  ScalarDimensionPair<float, 4>,
  ScalarDimensionPair<double, 1>,
  ScalarDimensionPair<double, 2>,
  ScalarDimensionPair<double, 3>,
  ScalarDimensionPair<double, 4>> UnsignedDistanceTypes;

typedef ::testing::Types<
  ScalarDimensionPair<float, 2>,
  ScalarDimensionPair<float, 3>,
  ScalarDimensionPair<float, 4>,
  ScalarDimensionPair<double, 2>,
  ScalarDimensionPair<double, 3>,
  ScalarDimensionPair<double, 4>> SignedDistanceTypes;

typedef ::testing::Types<
  ScalarDimensionPair<float, 2>,
  ScalarDimensionPair<double, 2>> AccuracyTypes;


// Fixtures.

template<typename T>
class UniformSpeedEikonalSolverTest : public ::testing::Test {
protected:
  virtual ~UniformSpeedEikonalSolverTest() {}
};


template<typename T>
class HighAccuracyUniformSpeedEikonalSolverTest : public ::testing::Test {
protected:
  virtual ~HighAccuracyUniformSpeedEikonalSolverTest() {}
};


template<typename T>
class VaryingSpeedEikonalSolverTest : public ::testing::Test {
protected:
  virtual ~VaryingSpeedEikonalSolverTest() {}
};


template<typename T>
class HighAccuracyVaryingSpeedEikonalSolverTest : public ::testing::Test {
protected:
  virtual ~HighAccuracyVaryingSpeedEikonalSolverTest() {}
};


template<typename T>
class UnsignedDistanceTest : public ::testing::Test {
protected:
  virtual ~UnsignedDistanceTest() {}
};


template<typename T>
class SignedDistanceTest : public ::testing::Test {
protected:
  virtual ~SignedDistanceTest() {}
};


template<typename T>
class UnsignedDistanceAccuracyTest : public ::testing::Test {
protected:
  virtual ~UnsignedDistanceAccuracyTest() {}
};


template<typename T>
class SignedDistanceAccuracyTest : public ::testing::Test {
protected:
  virtual ~SignedDistanceAccuracyTest() {}
};


// Associate types with fixtures.

TYPED_TEST_CASE(UniformSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(HighAccuracyUniformSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(VaryingSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(HighAccuracyVaryingSpeedEikonalSolverTest, EikonalSolverTypes);
TYPED_TEST_CASE(UnsignedDistanceTest, UnsignedDistanceTypes);
TYPED_TEST_CASE(SignedDistanceTest, SignedDistanceTypes);
TYPED_TEST_CASE(UnsignedDistanceAccuracyTest, AccuracyTypes);
TYPED_TEST_CASE(SignedDistanceAccuracyTest, AccuracyTypes);


// UniformSpeedEikonalSolverTest fixture.

TYPED_TEST(UniformSpeedEikonalSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 3>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element; // Invalid i'th element.
      auto const speed = ScalarType{1};

      auto expected_reason = stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [=]() {
          auto const eikonal_solver = EikonalSolverType(grid_spacing, speed);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(UniformSpeedEikonalSolverTest, InvalidSpeedThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = array<ScalarType, 3>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = invalid_speed; // Invalid speed!

    auto expected_reason = stringstream();
    expected_reason << "invalid speed: " << speed;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver = EikonalSolverType(grid_spacing, speed);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}


// HighAccuracyUniformSpeedEikonalSolverTest fixture.

TYPED_TEST(HighAccuracyUniformSpeedEikonalSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::HighAccuracyUniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 3>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element; // Invalid i'th element.
      auto const speed = ScalarType{1};

      auto expected_reason = stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [=]() {
          auto const eikonal_solver = EikonalSolverType(grid_spacing, speed);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(HighAccuracyUniformSpeedEikonalSolverTest, InvalidSpeedThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::HighAccuracyUniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = array<ScalarType, 3>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = invalid_speed; // Invalid speed!

    auto expected_reason = stringstream();
    expected_reason << "invalid speed: " << speed;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver = EikonalSolverType(grid_spacing, speed);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}


// VaryingSpeedEikonalSolverTest fixture.

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 3>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element; // Invalid i'th element.
      auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
      auto const speed_buffer =
        vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

      auto expected_reason = stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [&]() {
          auto const eikonal_solver =
            EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = array<ScalarType, 3>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

    // Invalid speed in the middle of the buffer.
    auto speed_buffer =
      vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});
    speed_buffer[speed_buffer.size() / 2] = invalid_speed;

    auto expected_reason = stringstream();
    expected_reason << "invalid speed: " << invalid_speed;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedBufferThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

  // Buffer size is not linear grid size!
  auto const speed_buffer =
    vector<ScalarType>(util::LinearSize(speed_grid_size) - 1, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "grid size " << util::ToString(speed_grid_size)
     << " does not match cell buffer size " << speed_buffer.size();

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const eikonal_solver =
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, InvalidSpeedGridSizeThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  for (auto i = size_t{0}; i < kDimension; ++i) {
    speed_grid_size[i] = 0; // Invalid i'th element.
    auto const speed_buffer =
      vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

    auto expected_reason = stringstream();
    expected_reason << "invalid grid size: "
                    << util::ToString(speed_grid_size);

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(VaryingSpeedEikonalSolverTest, IndexOutsideSpeedGridThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::VaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  // Speed grid smaller than distance grid!
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{9});
  auto const speed_buffer =
    vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  auto const boundary_distances = vector<ScalarType>(1, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "index outside speed grid";

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = thinks::fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}


// HighAccuracyVaryingSpeedEikonalSolverTest fixture.

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidGridSpacingThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_grid_spacing_elements = array<ScalarType, 3>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_grid_spacing_element : invalid_grid_spacing_elements)
  {
    for (auto i = size_t{0}; i < kDimension; ++i) {
      auto grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
      grid_spacing[i] = invalid_grid_spacing_element; // Invalid i'th element.
      auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
      auto const speed_buffer =
        vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

      auto expected_reason = stringstream();
      expected_reason << "invalid grid spacing: "
                      << util::ToString(grid_spacing);

      // Act.
      auto const ft = util::FunctionThrows<invalid_argument>(
        [&]() {
          auto const eikonal_solver =
            EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
        });

      // Assert.
      ASSERT_TRUE(ft.first);
      ASSERT_EQ(expected_reason.str(), ft.second);
    }
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidSpeedThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_speeds = array<ScalarType, 3>{{
    ScalarType{0},
    ScalarType{-1},
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_speed : invalid_speeds) {
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});

    // Invalid speed in the middle of the buffer.
    auto speed_buffer =
      vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});
    speed_buffer[speed_buffer.size() / 2] = invalid_speed;

    auto expected_reason = stringstream();
    expected_reason << "invalid speed: " << invalid_speed;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidSpeedBufferThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  // Buffer size is not the same as linear grid size!
  auto const speed_buffer =
    vector<ScalarType>(util::LinearSize(speed_grid_size) - 1, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "grid size " << util::ToString(speed_grid_size)
     << " does not match cell buffer size " << speed_buffer.size();

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const eikonal_solver =
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidSpeedGridSizeThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto speed_grid_size = util::FilledArray<kDimension>(size_t{10});
  for (auto i = size_t{0}; i < kDimension; ++i) {
    speed_grid_size[i] = 0; // Invalid i'th element.
    auto const speed_buffer =
      vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

    auto expected_reason = stringstream();
    expected_reason << "invalid grid size: " << util::ToString(speed_grid_size);

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const eikonal_solver =
          EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer);
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(HighAccuracyVaryingSpeedEikonalSolverTest, InvalidSpeedGridThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::HighAccuracyVaryingSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  // Speed grid smaller than distance grid!
  auto const speed_grid_size = util::FilledArray<kDimension>(size_t{9});
  auto const speed_buffer =
    vector<ScalarType>(util::LinearSize(speed_grid_size), ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  auto const boundary_distances = vector<ScalarType>(1, ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = thinks::fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed_grid_size, speed_buffer));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("index outside speed grid", ft.second);
}


// UnsignedDistance fixture.
// TODO - different dx.
TYPED_TEST(UnsignedDistanceTest, ZeroElementInGridSizeThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  for (auto i = size_t{0}; i < kDimension; ++i) {
    auto grid_size = util::FilledArray<kDimension>(size_t{10});
    grid_size[i] = 0; // Zero element in i'th position.
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = ScalarType{1};

    auto boundary_indices = vector<array<int32_t, kDimension>>();
    boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{0}));
    auto const boundary_distances = vector<ScalarType>(1, ScalarType{1});

    auto expected_reason = stringstream();
    expected_reason << "invalid grid size: " << util::ToString(grid_size);

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const unsigned_distance = thinks::fmm::UnsignedDistance(
          grid_size,
          boundary_indices,
          boundary_distances,
          EikonalSolverType(grid_spacing, speed));
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(UnsignedDistanceTest, EmptyBoundaryIndicesThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};
  auto const boundary_indices = vector<array<int32_t, kDimension>>{}; // Empty.
  auto const boundary_distances = vector<ScalarType>{}; // Empty.

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = thinks::fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("empty boundary condition", ft.second);
}

TYPED_TEST(UnsignedDistanceTest, FullGridBoundaryIndicesThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension> EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  // Add every cell in the grid!
  auto frozen_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  while (index_iter.has_next()) {
    frozen_indices.push_back(index_iter.index());
    index_iter.Next();
  }

  auto const frozen_distances =
    vector<ScalarType>(util::LinearSize(grid_size), ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = thinks::fmm::UnsignedDistance(
        grid_size,
        frozen_indices,
        frozen_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("whole grid is boundary", ft.second);
}

TYPED_TEST(UnsignedDistanceTest, DuplicateBoundaryIndicesThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension> EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  auto frozen_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  frozen_indices.push_back(index_iter.index());
  frozen_indices.push_back(index_iter.index()); // Same index!

  auto const frozen_distances = vector<ScalarType>(2, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "duplicate boundary index: "
                  << util::ToString(index_iter.index());

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const eikonal_solver =
        thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>(grid_spacing, speed);
      auto const unsigned_distance = thinks::fmm::UnsignedDistance(
        grid_size,
        frozen_indices,
        frozen_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(UnsignedDistanceTest, BoundaryIndexOutsideGridThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  auto frozen_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  frozen_indices.push_back(index_iter.index());
  frozen_indices.push_back(util::FilledArray<kDimension>(int32_t{-1})); // Outside!

  auto const frozen_distances = vector<ScalarType>(2, ScalarType{1});

  auto expected_reason = stringstream();
  expected_reason << "boundary index outside grid: "
                  << util::ToString(frozen_indices.back());

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = thinks::fmm::UnsignedDistance(
        grid_size,
        frozen_indices,
        frozen_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ(expected_reason.str(), ft.second);
}

TYPED_TEST(UnsignedDistanceTest, BoundaryIndicesAndDistancesSizeMismatchThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto index_iter = util::IndexIterator<kDimension>(grid_size);
  boundary_indices.push_back(index_iter.index());
  index_iter.Next();
  boundary_indices.push_back(index_iter.index());

  // Two indices, three distances.
  auto const boundary_distances = vector<ScalarType>(3, ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<invalid_argument>(
    [=]() {
      auto const unsigned_distance = thinks::fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("boundary indices/distances size mismatch", ft.second);
}

TYPED_TEST(UnsignedDistanceTest, InvalidBoundaryDistanceThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const invalid_boundary_distances = array<ScalarType, 3>{{
    ScalarType{-1},
    numeric_limits<ScalarType>::max(),
    numeric_limits<ScalarType>::quiet_NaN()
  }};

  for (auto const invalid_boundary_distance : invalid_boundary_distances) {
    auto const grid_size = util::FilledArray<kDimension>(size_t{10});
    auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
    auto const speed = ScalarType{1};

    auto boundary_indices = vector<array<int32_t, kDimension>>();
    auto index_iter = util::IndexIterator<kDimension>(grid_size);
    boundary_indices.push_back(index_iter.index());

    auto const boundary_distances =
      vector<ScalarType>(1, invalid_boundary_distance); // Invalid!

    auto expected_reason = stringstream();
    expected_reason << "invalid boundary distance: "
                    << invalid_boundary_distance;

    // Act.
    auto const ft = util::FunctionThrows<invalid_argument>(
      [=]() {
        auto const unsigned_distance = thinks::fmm::UnsignedDistance(
          grid_size,
          boundary_indices,
          boundary_distances,
          EikonalSolverType(grid_spacing, speed));
      });

    // Assert.
    ASSERT_TRUE(ft.first);
    ASSERT_EQ(expected_reason.str(), ft.second);
  }
}

TYPED_TEST(UnsignedDistanceTest, EikonalSolverFailThrows)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Eikonal solver cannot fail in one dimension.
  if (kDimension == 1) {
    return;
  }

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});
  auto const speed = ScalarType{1};

  // Create a scenario where solver has a very small value in one direction
  // and a very large value in another. Cannot resolve gradient for
  // this scenario.
  auto boundary_indices = vector<array<int32_t, kDimension>>();
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{0}));
  boundary_indices.push_back(util::FilledArray<kDimension>(int32_t{1}));
  auto boundary_distances = vector<ScalarType>();
  boundary_distances.push_back(ScalarType{1000});
  boundary_distances.push_back(ScalarType{1});

  // Act.
  auto const ft = util::FunctionThrows<runtime_error>(
    [=]() {
      auto const unsigned_distance = thinks::fmm::UnsignedDistance(
        grid_size,
        boundary_indices,
        boundary_distances,
        EikonalSolverType(grid_spacing, speed));
    });

  // Assert.
  ASSERT_TRUE(ft.first);
  ASSERT_EQ("negative discriminant", ft.second);
}

TYPED_TEST(UnsignedDistanceTest, DifferentUniformSpeed)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension> EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  auto frozen_indices = vector<array<int32_t, kDimension>>();
  frozen_indices.push_back(util::FilledArray<kDimension>(int32_t{5}));

  auto frozen_distances = vector<ScalarType>();
  frozen_distances.push_back(ScalarType{0});

  auto const speed1 = ScalarType{1};
  auto const speed2 = ScalarType{2};

  // Act.
  auto const unsigned_distance1 = thinks::fmm::UnsignedDistance(
    grid_size,
    frozen_indices,
    frozen_distances,
    EikonalSolverType(grid_spacing, speed1));

  auto const unsigned_distance2 = thinks::fmm::UnsignedDistance(
    grid_size,
    frozen_indices,
    frozen_distances,
    EikonalSolverType(grid_spacing, speed2));

  // Assert.
  // Check that the distance is halved when the speed is halved.
  // Note that the frozen distance is zero, which also passes this check.
  for (auto i = size_t{0}; i < unsigned_distance1.size(); ++i) {
    auto const d1 = unsigned_distance1[i];
    auto const d2 = unsigned_distance2[i];
    ASSERT_LE(fabs(speed1 * d1 - speed2 * d2), ScalarType(1e-3));
  }
}

TYPED_TEST(UnsignedDistanceTest, DISABLED_NonUniformGridSpacing)
{
  ASSERT_TRUE(false);
}


// UnsignedDistanceAccuracy fixture.

TYPED_TEST(UnsignedDistanceAccuracyTest, UniformSpeedHyperSphereStats)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{100});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.01));
  auto const uniform_speed = ScalarType{1};

  auto const sphere_center = util::FilledArray<kDimension>(ScalarType(0.5));
  auto const sphere_radius = ScalarType(0.25);

  auto const dist_stats =
    util::HyperSphereDistanceStats<ScalarType, kDimension>(
      sphere_center,
      sphere_radius,
      grid_size,
      grid_spacing,
      [](ScalarType const x) { return fabs(x); },
      [=](auto const frozen_indices,
          auto const frozen_distances) {
        return thinks::fmm::UnsignedDistance(
          grid_size,
          frozen_indices,
          frozen_distances,
          EikonalSolverType(grid_spacing, uniform_speed));
      });

  // Assert.
  typedef util::UnsignedDistanceAccuracyThresholds<ScalarType, kDimension>
    Thresholds;
  ASSERT_EQ(dist_stats.min_abs_error, Thresholds::min_abs_error_threshold());
  ASSERT_LT(dist_stats.max_abs_error, Thresholds::max_abs_error_threshold());
  ASSERT_LT(dist_stats.avg_abs_error, Thresholds::avg_abs_error_threshold());
  ASSERT_LT(dist_stats.std_dev_abs_error, Thresholds::std_dev_abs_error_threshold());
  //cerr << distance_stats << endl;
}

TYPED_TEST(UnsignedDistanceAccuracyTest, HighAccuracyUniformSpeedHyperSphereStats)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::HighAccuracyUniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  auto const grid_size = util::FilledArray<kDimension>(size_t{100});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.01));
  auto const uniform_speed = ScalarType{1};

  auto const sphere_center = util::FilledArray<kDimension>(ScalarType(0.5));
  auto const sphere_radius = ScalarType(0.25);

  auto const distance_stats =
    util::HyperSphereDistanceStats<ScalarType, kDimension>(
      sphere_center,
      sphere_radius,
      grid_size,
      grid_spacing,
      [](ScalarType const x) { return fabs(x); },
      [=](auto const frozen_indices,
          auto const frozen_distances) {
        return thinks::fmm::UnsignedDistance(
          grid_size,
          frozen_indices,
          frozen_distances,
          EikonalSolverType(grid_spacing, uniform_speed));
      });

  typedef util::UnsignedDistanceAccuracyThresholds<ScalarType, kDimension> Thresholds;
  ASSERT_EQ(distance_stats.min_abs_error, Thresholds::min_abs_error_threshold());
  ASSERT_LT(distance_stats.max_abs_error, Thresholds::max_abs_error_threshold());
  ASSERT_LT(distance_stats.avg_abs_error, Thresholds::avg_abs_error_threshold());
  ASSERT_LT(distance_stats.std_dev_abs_error, Thresholds::std_dev_abs_error_threshold());
  //cerr << "HIGH ACCURACY" << endl;
  //cerr << distance_stats << endl;
}

TYPED_TEST(UnsignedDistanceAccuracyTest, DISABLED_VaryingSpeed)
{
  ASSERT_TRUE(false);
}

TYPED_TEST(UnsignedDistanceAccuracyTest, DISABLED_HighAccuracyVaryingSpeed)
{
  ASSERT_TRUE(false);
}

TYPED_TEST(UnsignedDistanceTest, DISABLED_DistanceAccuracyOnHyperSphere)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension> EikonalSolverType;

  auto const grid_size = util::FilledArray<kDimension>(size_t{20});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.5));
  auto const speed = ScalarType{1};

  auto const center = util::FilledArray<kDimension>(ScalarType{5});
  auto const radius = ScalarType{3};
  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto boundary_distances = vector<ScalarType>();
  auto ground_truth = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    center, radius,
    grid_size, grid_spacing,
    [](ScalarType const d) { return fabs(d); },
    &boundary_indices, &boundary_distances,
    &ground_truth);

  auto const unsigned_distance = thinks::fmm::UnsignedDistance(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, speed));

  // TMP!!
  if (kDimension == 2) {
    // Negative values in shades of blue, positive values in shades of red.
    // Very large values as grey.
    auto const pixel_from_value = [](ScalarType const x) {
      if (x == numeric_limits<ScalarType>::max()) {
        return array<uint8_t, 3>{{128, 128, 128}};
      }
      return x < ScalarType(0) ?
        array<uint8_t, 3>{{
          uint8_t(0),
          uint8_t(0),
          util::Clamp<uint8_t>(
            ScalarType(0),
            ScalarType(numeric_limits<uint8_t>::max()),
            numeric_limits<uint8_t>::max() * fabs(x))}} :
        array<uint8_t, 3>{{
          util::Clamp<uint8_t>(
            ScalarType(0),
            ScalarType(numeric_limits<uint8_t>::max()),
            numeric_limits<uint8_t>::max() * x),
          uint8_t(0),
          uint8_t(0)}};
    };

    /*
    auto gt_values = vector<ScalarType>(
      grid_size[0] * grid_size[1], numeric_limits<ScalarType>::max());
    auto gt_grid = util::Grid<ScalarType, kDimension>(grid_size, gt_values.front());
    for (auto i = size_t{0}; i < frozen_indices.size(); ++i) {
      gt_grid.Cell(frozen_indices[i]) = frozen_distances[i];
    }
    */
    auto const gt_norm_values =
      util::SignedNormalized<ScalarType>(begin(ground_truth), end(ground_truth));
    auto ss_gt = stringstream();
    ss_gt << "gt_" << typeid(ScalarType).name() << ".ppm";
    thinks::ppm::writeRgbImage(
      ss_gt.str(),
      grid_size[0],
      grid_size[1],
      util::PixelsFromValues(
        begin(gt_norm_values),
        end(gt_norm_values),
        pixel_from_value));

    auto const usd_norm_values =
      util::SignedNormalized<ScalarType>(begin(unsigned_distance), end(unsigned_distance));
    auto ss_usd = stringstream();
    ss_usd << "unsigned_distance_" << typeid(ScalarType).name() << ".ppm";
    thinks::ppm::writeRgbImage(
      ss_usd.str(),
      grid_size[0],
      grid_size[1],
      util::PixelsFromValues(
        begin(usd_norm_values),
        end(usd_norm_values),
        pixel_from_value));
  }


  for (auto i = size_t{0}; i < unsigned_distance.size(); ++i) {
    if (kDimension == 2) {
      //cerr << "[" << i << "]: " << unsigned_distance[i] << " | " << ground_truth[i] << endl;
    }
    ASSERT_LE(fabs(unsigned_distance[i] - ground_truth[i]), 1e-3);
  }
}

TYPED_TEST(UnsignedDistanceAccuracyTest, GradientLength)
{
  using namespace std;

  typedef TypeParam::ScalarType ScalarType;
  static constexpr size_t kDimension = TypeParam::kDimension;
  typedef thinks::fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  auto const grid_size = util::FilledArray<kDimension>(size_t{100});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.01));
  auto const uniform_speed = ScalarType{1};

  auto const hyper_sphere_center = util::FilledArray<kDimension>(ScalarType(0.5));
  auto const hyper_sphere_radius = ScalarType(0.25);

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  auto boundary_distances = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    hyper_sphere_center,
    hyper_sphere_radius,
    grid_size,
    grid_spacing,
    [](ScalarType const x) { return fabs(x); },
    &boundary_indices,
    &boundary_distances);

#if 1
  auto const input_buffer =
    util::InputBuffer(grid_size, boundary_indices, boundary_distances);
  auto ss_input = stringstream();
  ss_input << "./grad_mag_input_" << typeid(ScalarType).name() << ".ppm";
  util::writeRgbImage(
    ss_input.str(),
    grid_size[0],
    grid_size[1],
    input_buffer);
#endif

  auto distance_buffer =
    thinks::fmm::UnsignedDistance(
      grid_size,
      boundary_indices,
      boundary_distances,
      EikonalSolverType(grid_spacing, uniform_speed));
  auto const distance_grid = util::Grid<ScalarType, kDimension>(
    grid_size, distance_buffer.front());

  auto const grad_buffer = util::DistanceGradients(distance_grid, grid_spacing);
  auto grad_mag_buffer = vector<ScalarType>(grad_buffer.size());
  for (auto i = size_t{0}; i < grad_buffer.size(); ++i) {
    grad_mag_buffer[i] = util::Magnitude(grad_buffer[i]);
  }
  auto const grad_mag_grid =
    util::Grid<ScalarType, kDimension>(grid_size, grad_mag_buffer.front());

  auto const grad_mag_error_buffer =
    util::GradientMagnitudeErrors(
      grad_mag_grid,
      boundary_indices);
#if 1
  auto ss_error = stringstream();
  ss_error << "./grad_mag_error_" << typeid(ScalarType).name() << ".ppm";
  util::writeRgbImage(
    ss_error.str(),
    grid_size[0],
    grid_size[1],
    grad_mag_error_buffer);
#endif

  auto const error_stats = util::ErrorStatistics(grad_mag_error_buffer);
  //cerr << "min: " << error_stats.min_abs_error << endl;
  //cerr << "max: " << error_stats.max_abs_error << endl;
  //cerr << "avg: " << error_stats.avg_abs_error << endl;

  //ASSERT_TRUE(false);
}

TYPED_TEST(UnsignedDistanceAccuracyTest, DISABLED_HighAccuracyGradientLength)
{
  ASSERT_TRUE(false);
}


// SignedDistance fixture.
// TODO


// SignedDistanceAccuracy fixture.
// TODO

#if 0

namespace detail {


template<typename T, typename R, typename U> inline
std::vector<R> TransformedVector(std::vector<T> const& v, U const unary_op)
{
  using namespace std;

  auto r = vector<R>(v.size());
  transform(begin(v), end(v), begin(r), unary_op);
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







template<typename T, std::size_t N> inline
std::vector<std::array<T, N>> DistanceGradients(
  std::vector<T>& distance_buffer,
  std::array<std::size_t, N> const& grid_size,
  std::array<T, N> const& grid_spacing)
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
      grid_spacing,
      Neighborhood<N>::offsets());
    valid_index = index_iter.Next();
  }

  return grad_buffer;
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
  auto const grid_spacing = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  HyperSphereFrozenCells<T, N>(
    center,
    radius,
    grid_size,
    grid_spacing,
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
    grid_spacing,
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
    grid_spacing);

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
  auto const grid_spacing = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  HyperSphereFrozenCells<T, N>(
    center,
    radius,
    grid_size,
    grid_spacing,
    &frozen_indices,
    &frozen_distances,
    &normals);

  auto const start_time = chrono::system_clock::now();
  auto distance_buffer = SignedDistance<T, N>(
    grid_size,
    grid_spacing,
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
    grid_spacing);

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
  auto const grid_spacing = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  auto distance_ground_truth_buffer = vector<T>();
  HyperSphereFrozenCells<T, N>(
    center,
    radius,
    grid_size,
    grid_spacing,
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
    grid_spacing,
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
  auto const grid_spacing = FilledArray<N>(T(0.01));
  auto const speed = T(1);

  auto frozen_indices = vector<array<int32_t, N>>();
  auto frozen_distances = vector<T>();
  auto normals = vector<array<T, N>>();
  auto distance_ground_truth_buffer = vector<T>();
  HyperSphereFrozenCells<T, N>(
    center,
    radius,
    grid_size,
    grid_spacing,
    &frozen_indices,
    &frozen_distances,
    &normals,
    &distance_ground_truth_buffer);

  auto const start_time = chrono::system_clock::now();
  auto distance_buffer = SignedDistance<T, N>(
    grid_size,
    grid_spacing,
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



#endif

} // namespace

namespace std {

template<typename T, size_t N>
ostream& operator<<(ostream& os, util::DistanceStats<T, N> const& dist_stats)
{
  os << "Distance value stats <" << typeid(T).name() << ", " << N << ">:" << endl
    << "min abs error: " << dist_stats.min_abs_error << endl
    << "max abs error: " << dist_stats.max_abs_error << endl
    << "avg abs error: " << dist_stats.avg_abs_error << endl
    << "std_dev abs error: " << dist_stats.std_dev_abs_error << endl
    << "duration: " << dist_stats.duration_in_s << " [s]" << endl;
  return os;
}

} // namespace std
