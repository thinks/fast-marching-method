#include <iostream>

#include <thinks/ppm.hpp>
#include "../include/thinks/fast_marching_method/fast_marching_method.hpp"
#include "../test/util.hpp"

namespace {

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

template<typename T, typename InIter> inline
std::vector<T> MaxAbsNormalized(InIter const in_begin, InIter const in_end)
{
  using namespace std;

  auto max_abs_value = T{0};
  for (auto iter = in_begin; iter != in_end; ++iter) {
    auto const value = *iter;
    if (value == numeric_limits<T>::max()) {
      continue;
    }

    max_abs_value = max(max_abs_value, abs(value));
  }
  auto r = vector<T>{};
  transform(
    in_begin,
    in_end,
    back_inserter(r),
    [=](auto const value) {
      if (value == numeric_limits<T>::max()) {
        return value;
      }
      return value / max_abs_value;
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

  //auto const norm_pixel_data =
  //  SignedNormalized<T>(begin(pixel_data), end(pixel_data));
  auto const norm_pixel_data =
    MaxAbsNormalized<T>(begin(pixel_data), end(pixel_data));
  thinks::ppm::writeRgbImage(
    filename,
    width,
    height,
    PixelsFromValues(
      begin(norm_pixel_data),
      end(norm_pixel_data),
      pixel_from_value));
}

template<typename T>
void writeGreyImage(
  std::string const& filename,
  std::size_t width,
  std::size_t height,
  std::vector<T> const& pixel_data)
{
  using namespace std;

  auto const pixel_from_value = [](T const x) {
    auto const intensity = Clamp<uint8_t>(
      T{0},
      T(numeric_limits<uint8_t>::max()),
      numeric_limits<uint8_t>::max() * fabs(x));
    return array<uint8_t, 3>{{intensity, intensity, intensity}};
  };

  auto const max_abs_normalized_pixel_data =
    MaxAbsNormalized<T>(begin(pixel_data), end(pixel_data));
  thinks::ppm::writeRgbImage(
    filename,
    width,
    height,
    PixelsFromValues(
      begin(max_abs_normalized_pixel_data),
      end(max_abs_normalized_pixel_data),
      pixel_from_value));
}

// Returns a grid buffer where boundary distances have been set. Cells that
// are not on the boundary have the value std::numeric_limits<T>::max().
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
      vector<T>(util::LinearSize(grid_size), numeric_limits<T>::max());
  auto input_grid = util::Grid<T, N>(grid_size, input_buffer.front());
  for (auto i = size_t{0}; i < boundary_indices.size(); ++i) {
    auto const boundary_index = boundary_indices[i];
    /*if (!util::Inside(boundary_index, grid_size)) {
      throw runtime_error("index outside grid");
    }*/

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
    [=](auto const v, auto const gt) { return (v - gt) / min_grid_spacing; });
  return r;
}


struct ErrorStats
{
  double min_abs_error;
  double max_abs_error;
  double avg_abs_error;
  double std_dev_abs_error;
};

std::ostream& operator<<(std::ostream& os, ErrorStats const& es)
{
  using namespace std;

  os << "min_abs_error: " << es.min_abs_error << endl;
  os << "max_abs_error: " << es.max_abs_error << endl;
  os << "avg_abs_error: " << es.avg_abs_error << endl;
  os << "std_dev_abs_error: " << es.std_dev_abs_error << endl;
  return os;
}


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
  util::Grid<T, N> const& grid,
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
  util::Grid<T, N> const& distance_grid,
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


//! Compute gradient magnitide errors on the form:
//! error = grad_mag - 1
template<typename T, std::size_t N>
std::vector<T> GradientMagnitudeErrors(
  util::Grid<T, N> const& grad_mag_grid,
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




#if 0
void TimingTest2D()
{
  using namespace std;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<double, 2>
    EikonalSolverType;
  typedef fmm::HighAccuracyUniformSpeedEikonalSolver<double, 2>
    HighAccuracyEikonalSolverType;
  typedef fmm::DistanceSolver<double, 2>
    DistanceSolverType;

  // Arrange.
  auto const grid_spacing = util::FilledArray<2>(double{1});
  auto const uniform_speed = double{1};

  // 2D
  auto eikonal_times = vector<double>();
  auto high_accuracy_eikonal_times = vector<double>();
  auto distance_times = vector<double>();
  eikonal_times.reserve(512);
  high_accuracy_eikonal_times.reserve(512);
  distance_times.reserve(512);
  for (auto i = size_t{64}; i <= 512; i += 2) {
    auto const grid_size = util::FilledArray<2>(i);

    // Simple point boundary for regular fast marching.
    auto boundary_indices = vector<array<int32_t, 2>>(
      size_t{1}, util::FilledArray<2>(int32_t(i / 2)));
    auto boundary_distances = vector<double>(size_t{1}, double{0});

    // Compute exact distances in vertex neighborhood for high accuracy
    // fast marching.
    auto high_accuracy_boundary_indices = vector<array<int32_t, 2>>();
    high_accuracy_boundary_indices.push_back(
      util::FilledArray<2>(int32_t(i / 2))); // Center.
    auto const vtx_neighbor_offsets = util::VertexNeighborOffsets<2>();
    for (auto const& vtx_neighbor_offset : vtx_neighbor_offsets) {
      auto index = boundary_indices[0];
      for (auto i = size_t{0}; i < 2; ++i) {
        index[i] += vtx_neighbor_offset[i];
      }
      high_accuracy_boundary_indices.push_back(index);
    }
    auto high_accuracy_boundary_distances = vector<double>();
    high_accuracy_boundary_distances.push_back(double{0}); // Center.
    auto center_position = util::FilledArray<2>(double{0});
    for (auto i = size_t{0}; i < 2; ++i) {
      center_position[i] =
        (high_accuracy_boundary_indices[0][i] + double(0.5)) * grid_spacing[i];
    }
    for (auto j = size_t{1}; j < high_accuracy_boundary_indices.size(); ++j) {
      auto const& index = high_accuracy_boundary_indices[j];
      auto position = util::FilledArray<2>(double(0));
      for (auto i = size_t{0}; i < 2; ++i) {
        position[i] = (index[i] + double(0.5)) * grid_spacing[i];
      }
      auto delta = util::FilledArray<2>(double(0));
      for (auto i = size_t{0}; i < 2; ++i) {
        delta[i] = center_position[i] - position[i];
      }
      high_accuracy_boundary_distances.push_back(util::Magnitude(delta));
    }

    {
    auto const start_time = chrono::system_clock::now();
    auto eikonal = fmm::UnsignedArrivalTime(
      grid_size,
      boundary_indices,
      boundary_distances,
      EikonalSolverType(grid_spacing, uniform_speed));
    auto const end_time = chrono::system_clock::now();
    eikonal_times.push_back(chrono::duration<double>(end_time - start_time).count());
    }

    {
    auto const start_time = chrono::system_clock::now();
    auto high_accuracy_eikonal = fmm::UnsignedArrivalTime(
      grid_size,
      high_accuracy_boundary_indices,
      high_accuracy_boundary_distances,
      HighAccuracyEikonalSolverType(grid_spacing, uniform_speed));
    auto const end_time = chrono::system_clock::now();
    high_accuracy_eikonal_times.push_back(chrono::duration<double>(end_time - start_time).count());
    }

    {
    auto const start_time = chrono::system_clock::now();
    auto distance = fmm::UnsignedArrivalTime(
      grid_size,
      boundary_indices,
      boundary_distances,
      DistanceSolverType(grid_spacing[0]));
    auto const end_time = chrono::system_clock::now();
    distance_times.push_back(chrono::duration<double>(end_time - start_time).count());
    }

    cerr << util::ToString(grid_size) << endl;
  }

  auto ofs = ofstream();
  ofs.open("./timing.txt", ofstream::trunc);
  ofs << "Eikonal" << endl;
  for (auto const t : eikonal_times) {
    ofs << t << endl;
  }
  ofs << "High accuracy Eikonal" << endl;
  for (auto const t : high_accuracy_eikonal_times) {
    ofs << t << endl;
  }
  ofs << "Distance" << endl;
  for (auto const t : distance_times) {
    ofs << t << endl;
  }
  ofs.close();
}
#endif

#if 0
void OverlappingCircles()
{
  using namespace std;

  typedef float ScalarType;
  static constexpr size_t kDimension = 2;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{50});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.02));
  auto const uniform_speed = ScalarType{1};

  auto sphere_center1 = util::FilledArray<kDimension>(ScalarType{0.5});
  sphere_center1[0] = ScalarType(0.4);
  auto const sphere_radius1 = ScalarType{0.25};
  auto sphere_center2 = util::FilledArray<kDimension>(ScalarType{0.5});
  sphere_center2[0] = ScalarType(0.6);
  auto const sphere_radius2 = ScalarType{0.25};

  auto sphere_boundary_indices1 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances1 = vector<ScalarType>();
  util::BoxBoundaryCells(
    util::FilledArray<kDimension>(int32_t{5}),
    util::FilledArray<kDimension>(size_t{20}),
    grid_size,
    &sphere_boundary_indices1,
    &sphere_boundary_distances1);
  /*util::HyperSphereBoundaryCells(
    sphere_center1,
    sphere_radius1,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return fabs(d); },
    &sphere_boundary_indices1,
    &sphere_boundary_distances1);*/
  auto sphere_boundary_indices2 = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_distances2 = vector<ScalarType>();
  util::BoxBoundaryCells(
    util::FilledArray<kDimension>(int32_t{15}),
    util::FilledArray<kDimension>(size_t{20}),
    grid_size,
    &sphere_boundary_indices2,
    &sphere_boundary_distances2);
  /*util::HyperSphereBoundaryCells(
    sphere_center2,
    sphere_radius2,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return fabs(d); },
    &sphere_boundary_indices2,
    &sphere_boundary_distances2);*/

  auto boundary_indices = sphere_boundary_indices1;
  auto boundary_distances = sphere_boundary_distances1;
  for (auto i = size_t{0}; i < sphere_boundary_indices2.size(); ++i) {
    auto const index = sphere_boundary_indices2[i];
    if (find(begin(boundary_indices), end(boundary_indices), index) ==
        end(boundary_indices)) {
      boundary_indices.push_back(index);
      boundary_distances.push_back(sphere_boundary_distances2[i]);
    }
  }

  auto input = InputBuffer(grid_size, boundary_indices, boundary_distances);
  writeRgbImage("./overlap_input.ppm", grid_size[0], grid_size[1], input);

  // Act.
  auto const unsigned_distance = fmm::UnsignedArrivalTime(
    grid_size,
    boundary_indices,
    boundary_distances,
    EikonalSolverType(grid_spacing, uniform_speed));
  writeRgbImage("./overlap_result.ppm", grid_size[0], grid_size[1], unsigned_distance);
}
#endif

#if 0
void CircleError()
{
  using namespace std;

  typedef float ScalarType;
  static constexpr size_t kDimension = 2;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{50});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.02));
  auto const uniform_speed = ScalarType{1};

  auto const sphere_center = util::FilledArray<kDimension>(ScalarType{0.5});
  auto const sphere_radius = ScalarType{0.25};

  auto unsigned_sphere_boundary_indices = vector<array<int32_t, kDimension>>();
  auto unsigned_sphere_boundary_times = vector<ScalarType>();
  auto unsigned_gt = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    sphere_center,
    sphere_radius,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return fabs(d); },
    &unsigned_sphere_boundary_indices,
    &unsigned_sphere_boundary_times,
    &unsigned_gt);
  auto signed_sphere_boundary_indices = vector<array<int32_t, kDimension>>();
  auto signed_sphere_boundary_times = vector<ScalarType>();
  auto signed_gt = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    sphere_center,
    sphere_radius,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return d; },
    &signed_sphere_boundary_indices,
    &signed_sphere_boundary_times,
    &signed_gt);

  auto unsigned_input = InputBuffer(grid_size, unsigned_sphere_boundary_indices, unsigned_sphere_boundary_times);
  writeRgbImage("./circle_unsigned_input.ppm", grid_size[0], grid_size[1], unsigned_input);
  auto signed_input = InputBuffer(grid_size, signed_sphere_boundary_indices, signed_sphere_boundary_times);
  writeRgbImage("./circle_signed_input.ppm", grid_size[0], grid_size[1], signed_input);

  // Act.
  auto const unsigned_times = fmm::UnsignedArrivalTime(
    grid_size,
    unsigned_sphere_boundary_indices,
    unsigned_sphere_boundary_times,
    EikonalSolverType(grid_spacing, uniform_speed));
  writeRgbImage("./circle_unsigned_result.ppm", grid_size[0], grid_size[1], unsigned_times);

  auto const signed_times = fmm::SignedArrivalTime(
    grid_size,
    signed_sphere_boundary_indices,
    signed_sphere_boundary_times,
    EikonalSolverType(grid_spacing, uniform_speed));
  writeRgbImage("./circle_signed_result.ppm", grid_size[0], grid_size[1], signed_times);

  auto unsigned_errors = ErrorBuffer(unsigned_times, unsigned_gt, grid_spacing);
  writeGreyImage("./circle_unsigned_error.ppm", grid_size[0], grid_size[1], unsigned_errors);
  auto signed_errors = ErrorBuffer(signed_times, signed_gt, grid_spacing);
  writeGreyImage("./circle_signed_error.ppm", grid_size[0], grid_size[1], signed_errors);

  auto unsigned_error_stats = ErrorStatistics(unsigned_errors);
  cout << "unsigned_error_stats" << endl;
  cout << unsigned_error_stats << endl;
  auto signed_error_stats = ErrorStatistics(signed_errors);
  cout << "signed_error_stats" << endl;
  cout << signed_error_stats << endl;
}
#endif

void CheckerBoard()
{
  using namespace std;

  typedef float ScalarType;
  static constexpr auto kDimension = 2;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  // Arrange.
  auto const grid_size = util::FilledArray<kDimension>(size_t{10});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType{1});

  auto boundary_indices = vector<array<int32_t, kDimension>>();
  {
    auto index_iter = util::IndexIterator<kDimension>(grid_size);
    auto even = [](auto const i) { return i % 2 == 0; };
    while (index_iter.has_next()) {
      auto const index = index_iter.index();
      if (all_of(begin(index), end(index), even) ||
          none_of(begin(index), end(index), even)) {
         boundary_indices.push_back(index);
      }
      index_iter.Next();
    }
  }

  auto boundary_times = vector<ScalarType>(
    boundary_indices.size(), ScalarType{0});

  auto const speed = ScalarType{1};

  // Act.
  auto signed_times = fmm::SignedArrivalTime(
    grid_size,
    boundary_indices,
    boundary_times,
    EikonalSolverType(grid_spacing, speed));

  writeRgbImage("./checkerboard_result.ppm", grid_size[0], grid_size[1], signed_times);


#if 0
  // Assert.
  auto time_grid = util::Grid<ScalarType, kDimension>(
    grid_size, signed_times.front());
  {
    auto index_iter = util::IndexIterator<kDimension>(grid_size);
    auto const odd = [](auto const i) { return i % 2 == 1; };
    while (index_iter.has_next()) {
      auto const index = index_iter.index();

      auto is_boundary = false;
      for (auto i = size_t{0}; i < kDimension; ++i) {
        if (index[i] == 0 || index[i] == grid_size[i] - 1) {
          is_boundary = true;
          break;
        }
      }

      auto const time = time_grid.Cell(index);
      if (!is_boundary &&
          (all_of(begin(index), end(index), odd) ||
           none_of(begin(index), end(index), odd))) {
        ASSERT_LE(time, ScalarType{0});
      }
      else {
        ASSERT_GE(time, ScalarType{0});
      }
      index_iter.Next();
    }
  }
#endif
}


void InputConcept()
{
  using namespace std;

  typedef float ScalarType;
  static constexpr size_t kDimension = 2;
  namespace fmm = thinks::fast_marching_method;
  typedef fmm::UniformSpeedEikonalSolver<ScalarType, kDimension>
    EikonalSolverType;

  auto const grid_size = util::FilledArray<kDimension>(size_t{1000});
  auto const grid_spacing = util::FilledArray<kDimension>(ScalarType(0.001));
  auto const uniform_speed = ScalarType{1};

  auto const sphere_center = util::FilledArray<kDimension>(ScalarType{0.5});
  auto const sphere_radius = ScalarType{0.25f * 1.12f};

  auto const dilation_pass_count = size_t{0};
  auto sphere_boundary_indices = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_times = vector<ScalarType>();
  auto gt = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    sphere_center,
    sphere_radius,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return d; },
    dilation_pass_count,
    &sphere_boundary_indices,
    &sphere_boundary_times,
    &gt);

  auto const arrival_times = fmm::SignedArrivalTime(
    grid_size,
    sphere_boundary_indices,
    sphere_boundary_times,
    EikonalSolverType(grid_spacing, uniform_speed));
  writeRgbImage("./input_figure.ppm", grid_size[0], grid_size[1], arrival_times);
}

void InputExample()
{
  using namespace std;
  namespace fmm = thinks::fast_marching_method;

  /*
   * Generated using the following...
   *
  auto sphere_center = {{0.5f, 0.5f}};
  auto sphere_radius = 0.28f;
  auto const dilation_pass_count = size_t{0};
  auto sphere_boundary_indices = vector<array<int32_t, kDimension>>();
  auto sphere_boundary_times = vector<ScalarType>();
  auto gt = vector<ScalarType>();
  util::HyperSphereBoundaryCells(
    sphere_center,
    sphere_radius,
    grid_size,
    grid_spacing,
    [](ScalarType const d) { return d; },
    dilation_pass_count,
    &sphere_boundary_indices,
    &sphere_boundary_times,
    &gt);
  */
  auto circle_boundary_indices = vector<array<int32_t, 2>>{
    {{5, 3}}, {{6, 3}}, {{7, 3}}, {{8, 3}}, {{9, 3}}, {{10, 3}}, {{4, 4}},
    {{5, 4}}, {{10, 4}}, {{11, 4}}, {{3, 5}}, {{4, 5}}, {{11, 5}}, {{12, 5}},
    {{3, 6}}, {{12, 6}}, {{3, 7}}, {{12, 7}}, {{3, 8}}, {{12, 8}}, {{3, 9}},
    {{12, 9}}, {{3, 10}}, {{4, 10}}, {{11, 10}}, {{12, 10}}, {{4, 11}},
    {{5, 11}}, {{10, 11}}, {{11, 11}}, {{5, 12}}, {{6, 12}}, {{7, 12}},
    {{8, 12}}, {{9, 12}}, {{10, 12}},
  };
  auto circle_boundary_times = vector<float>{
    0.0417385f, 0.0164635f, 0.0029808f, 0.0029808f, 0.0164635f, 0.0417385f,
    0.0293592f, -0.0111773f, -0.0111773f, 0.0293592f, 0.0417385f, -0.0111773f,
    -0.0111773f, 0.0417385f, 0.0164635f, 0.0164635f, 0.0029808f, 0.0029808f,
    0.0029808f, 0.0029808f, 0.0164635f, 0.0164635f, 0.0417385f, -0.0111773f,
    -0.0111773f, 0.0417385f, 0.0293592f, -0.0111773f, -0.0111773f, 0.0293592f,
    0.0417385f, 0.0164635f, 0.0029808f, 0.0029808f, 0.0164635f, 0.0417385f
  };

  auto grid_size = array<size_t, 2>{{16, 16}};
  auto grid_spacing = array<float, 2>{{1.f/16, 1.f/16}};
  auto uniform_speed = 1.f;

  auto arrival_times = fmm::SignedArrivalTime(
    grid_size,
    circle_boundary_indices,
    circle_boundary_times,
    fmm::UniformSpeedEikonalSolver<float, 2>(grid_spacing, uniform_speed));
  writeRgbImage("./input_figure_values.ppm", grid_size[0], grid_size[1], arrival_times);

  /*
  for (auto i = 0; i < sphere_boundary_indices.size(); ++i) {
    auto const index = sphere_boundary_indices[i];
    auto const time = sphere_boundary_times[i];
    cout << "[" << index[0] << ", " << index[1] << "] = " << time << endl;
  }
  */

#if 0
  auto input_buffer = InputBuffer(grid_size, circle_boundary_indices, circle_boundary_times);
  writeRgbImage("./input_figure_values_input.ppm", grid_size[0], grid_size[1], input_buffer);
#endif
}

void EikonalSolversPointSource()
{
  using namespace std;
  using namespace util;
  namespace fmm = thinks::fast_marching_method;

  typedef double S;

  auto add = [](array<int32_t, 2> const& u, array<int32_t, 2> const& v)
  {
    return array<int32_t, 2>{{u[0] + v[0], u[1] + v[1]}};
  };

  auto const grid_size = array<size_t, 2>{{
    512,
    512
  }};
  auto const grid_spacing = array<S, 2>{{
    S{1} / grid_size[0],
    S{1} / grid_size[1]
  }};

  auto sphere_center = array<S, 2>{{S(0.5), S(0.5)}};
  auto sphere_radius = S(0.01);

  auto const dilation_pass_count = size_t{0};
  auto boundary_indices = vector<array<int32_t, 2>>();
  auto boundary_times = vector<S>();
  auto gt = vector<S>();
  util::HyperSphereBoundaryCells(
    sphere_center,
    sphere_radius,
    grid_size,
    grid_spacing,
    [](S const d) { return d; },
    dilation_pass_count,
    &boundary_indices,
    &boundary_times,
    &gt);
  auto const ha_dilation_pass_count = size_t{1};
  auto ha_boundary_indices = vector<array<int32_t, 2>>();
  auto ha_boundary_times = vector<S>();
  auto ha_gt = vector<S>();
  util::HyperSphereBoundaryCells(
    sphere_center,
    sphere_radius,
    grid_size,
    grid_spacing,
    [](S const d) { return d; },
    ha_dilation_pass_count,
    &ha_boundary_indices,
    &ha_boundary_times,
    &ha_gt);

  /*
  auto point_boundary_indices = vector<array<int32_t, 2>>{
    {{256, 256}}
  };
  auto point_boundary_times = vector<float>{0.f};
  auto ha_point_boundary_indices = vector<array<int32_t, 2>>{
    point_boundary_indices[0],
    add(point_boundary_indices[0], {{-1, 0}}),
    add(point_boundary_indices[0], {{+1, 0}}),
    add(point_boundary_indices[0], {{0, -1}}),
    add(point_boundary_indices[0], {{0, +1}}),
    add(point_boundary_indices[0], {{-1, -1}}),
    add(point_boundary_indices[0], {{-1, +1}}),
    add(point_boundary_indices[0], {{+1, -1}}),
    add(point_boundary_indices[0], {{+1, +1}}),
  };
  auto ha_point_boundary_times = vector<float>{
    0.f,
    Distance(
      CellCenter(ha_point_boundary_indices[0], grid_spacing),
      CellCenter(ha_point_boundary_indices[1], grid_spacing)),
    Distance(
      CellCenter(ha_point_boundary_indices[0], grid_spacing),
      CellCenter(ha_point_boundary_indices[2], grid_spacing)),
    Distance(
      CellCenter(ha_point_boundary_indices[0], grid_spacing),
      CellCenter(ha_point_boundary_indices[3], grid_spacing)),
    Distance(
      CellCenter(ha_point_boundary_indices[0], grid_spacing),
      CellCenter(ha_point_boundary_indices[4], grid_spacing)),
    Distance(
      CellCenter(ha_point_boundary_indices[0], grid_spacing),
      CellCenter(ha_point_boundary_indices[5], grid_spacing)),
    Distance(
      CellCenter(ha_point_boundary_indices[0], grid_spacing),
      CellCenter(ha_point_boundary_indices[6], grid_spacing)),
    Distance(
      CellCenter(ha_point_boundary_indices[0], grid_spacing),
      CellCenter(ha_point_boundary_indices[7], grid_spacing)),
    Distance(
      CellCenter(ha_point_boundary_indices[0], grid_spacing),
      CellCenter(ha_point_boundary_indices[8], grid_spacing)),
  };
  */

  // Uniform speed.
  auto const uniform_speed = S(1);
  {
    auto arrival_times = fmm::SignedArrivalTime(
      grid_size,
      boundary_indices,
      boundary_times,
      fmm::UniformSpeedEikonalSolver<S, 2>(
        grid_spacing,
        uniform_speed));
    writeRgbImage(
      "./uniform_speed_input.ppm",
      grid_size[0],
      grid_size[1],
      InputBuffer(grid_size, boundary_indices, boundary_times));
    writeRgbImage(
      "./uniform_speed.ppm",
      grid_size[0],
      grid_size[1],
      arrival_times);
  }

  {
    auto arrival_times = fmm::SignedArrivalTime(
      grid_size,
      ha_boundary_indices,
      ha_boundary_times,
      fmm::HighAccuracyUniformSpeedEikonalSolver<S, 2>(
        grid_spacing,
        uniform_speed));
    writeRgbImage(
      "./ha_uniform_speed_input.ppm",
      grid_size[0],
      grid_size[1],
      InputBuffer(grid_size, ha_boundary_indices, ha_boundary_times));
    writeRgbImage(
      "./ha_uniform_speed.ppm",
      grid_size[0],
      grid_size[1],
      arrival_times);
  }

  // Varying speed.
  auto const speed_grid_size = grid_size;
  auto speed_grid_buffer = vector<S>(
    LinearSize(grid_size),
    S{1});
  auto speed_grid = Grid<S, 2>(speed_grid_size, speed_grid_buffer.front());
  auto speed_grid_iter = IndexIterator<2>(speed_grid.size());
  while (speed_grid_iter.has_next()) {
    auto const index = speed_grid_iter.index();
    if (index[0] < static_cast<int32_t>(speed_grid.size()[0] / 4)) {
      speed_grid.Cell(index) = 1.5f;
    }
    speed_grid_iter.Next();
  }

  {
    auto arrival_times = fmm::SignedArrivalTime(
      grid_size,
      boundary_indices,
      boundary_times,
      fmm::VaryingSpeedEikonalSolver<S, 2>(
        grid_spacing,
        speed_grid_size,
        speed_grid_buffer));
    writeRgbImage(
      "./varying_speed_input.ppm",
      grid_size[0],
      grid_size[1],
      InputBuffer(grid_size, boundary_indices, boundary_times));
    writeRgbImage(
      "./varying_speed.ppm",
      grid_size[0],
      grid_size[1],
      arrival_times);
  }

  {
    auto arrival_times = fmm::SignedArrivalTime(
      grid_size,
      ha_boundary_indices,
      ha_boundary_times,
      fmm::HighAccuracyVaryingSpeedEikonalSolver<S, 2>(
        grid_spacing,
        speed_grid_size,
        speed_grid_buffer));
    writeRgbImage(
      "./ha_varying_speed_input.ppm",
      grid_size[0],
      grid_size[1],
      InputBuffer(grid_size, ha_boundary_indices, ha_boundary_times));
    writeRgbImage(
      "./ha_varying_speed.ppm",
      grid_size[0],
      grid_size[1],
      arrival_times);
  }
}

} // namespace


int main(int argc, char* argv[])
{
  using namespace std;

  //OverlappingCircles();
  //CircleError();
  //CheckerBoard();
  try {
    InputConcept();
    InputExample();
    //EikonalSolversPointSource();

    cin.get(); // Keep console open...
    return 0;
  }
  catch (std::exception& ex) {
    cerr << ex.what() << endl;
    cin.get(); // Close console.
    return 1; 
  }
}
