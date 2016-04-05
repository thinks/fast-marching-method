#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <thinks/testFastMarchingMethod.hpp>
#include <thinks/ppm.hpp>


namespace {

template<typename R, typename T> inline
R clamp(T const low, T const high, T const value)
{
  using namespace std;

  return static_cast<R>(min<T>(high, max<T>(low, value)));
}

template<typename R, typename U, typename InIter> inline
std::vector<R> transformedAsVector(
  InIter const begin, InIter const end, U const unary_op)
{
  using namespace std;

  auto r = vector<R>();
  transform(begin, end, back_inserter(r), unary_op);
  return r;
}

template<typename T, typename InIter> inline
std::vector<T> signedNormalized(InIter const in_begin, InIter const in_end)
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

  return transformedAsVector<T>(
    in_begin,
    in_end,
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
}

template<typename InIter, typename C> inline
std::vector<std::uint8_t> pixelsFromValues(
  InIter const in_begin,
  InIter const in_end,
  C const pixel_from_value)
{
  using namespace std;

  auto pixels = vector<uint8_t>();
  for (auto iter = in_begin; iter != in_end; ++iter) {
    auto const pixel = pixel_from_value(*iter);
    pixels.insert(end(pixels), begin(pixel), end(pixel));
  }
  return pixels;
}

template<typename T>
void writeGradMagImages(
  thinks::fmm::test::GradientMagnitudeStats<T, 2> const& grad_mag_stats,
  std::string const& prefix)
{
  using namespace std;
  using namespace thinks;

  // Negative values in shades of blue, positive values in shades of red.
  // Very large values as grey.
  auto const pixel_from_value = [](T const x) {
    if (x == numeric_limits<T>::max()) {
      return array<uint8_t, 3>{{128, 128, 128}};
    }
    return x < T(0) ?
      array<uint8_t, 3>{{
        uint8_t(0),
        uint8_t(0),
        clamp<uint8_t>(
          T(0),
          T(numeric_limits<uint8_t>::max()),
          numeric_limits<uint8_t>::max() * fabs(x))}} :
      array<uint8_t, 3>{{
        clamp<uint8_t>(
          T(0),
          T(numeric_limits<uint8_t>::max()),
          numeric_limits<uint8_t>::max() * x),
        uint8_t(0),
        uint8_t(0)}};
  };

  auto const width = grad_mag_stats.grid_size[0];
  auto const height = grad_mag_stats.grid_size[1];

  stringstream ss_input;
  ss_input << prefix << "_input_" << typeid(T).name() << ".ppm";
  auto const normalized_input = signedNormalized<T>(
    begin(grad_mag_stats.input_buffer),
    end(grad_mag_stats.input_buffer));
  ppm::writeRgbImage(
    ss_input.str(),
    width,
    height,
    pixelsFromValues(
      begin(normalized_input),
      end(normalized_input),
      pixel_from_value));

  stringstream ss_distance;
  ss_distance << prefix << "_distance_" << typeid(T).name() << ".ppm";
  auto const normalized_distance = signedNormalized<T>(
    begin(grad_mag_stats.distance_buffer),
    end(grad_mag_stats.distance_buffer));
  ppm::writeRgbImage(
    ss_distance.str(),
    width,
    height,
    pixelsFromValues(
      begin(normalized_distance),
      end(normalized_distance),
      pixel_from_value));

#if 0
  stringstream ss_grad_mag;
  ss_grad_mag << prefix << "_" << typeid(T).name() << ".ppm";
  auto const grad_mag = transformedAsVector<T>(
    begin(grad_mag_stats.grad_buffer),
    end(grad_mag_stats.grad_buffer),
    [](auto const v) { return sqrt(v[0] * v[0] + v[1] * v[1]); } ),
  ppm::writeRgbImage(
    ss_grad_mag.str(),
    width,
    height,
    pixelsFromValues(
      begin(grad_mag),
      end(grad_mag),
      [](T const& v) {
        auto const
        return array<uint8_t, 3>{{}};
          clamp<Pixel8::ChannelType>(
            T(0),
            T(numeric_limits<Pixel8::ChannelType>::max()),
            numeric_limits<Pixel8::ChannelType>::max() * v));
      }));
#endif

  stringstream ss_error;
  ss_error << prefix << "_error_" << typeid(T).name() << ".ppm";
  auto const normalized_error = signedNormalized<T>(
    begin(grad_mag_stats.error_buffer),
    end(grad_mag_stats.error_buffer));
  ppm::writeRgbImage(
    ss_error.str(),
    width,
    height,
    pixelsFromValues(
      begin(normalized_error),
      end(normalized_error),
      pixel_from_value));
}

template<typename T>
void writeDistStatImages(
  thinks::fmm::test::DistanceValueStats<T, 2> const& dist_stats,
  std::string const& prefix)
{
  using namespace std;
  using namespace thinks;

  // Negative values in shades of blue, positive values in shades of red.
  // Very large values as grey.
  auto const pixel_from_value = [](T const x) {
    if (x == numeric_limits<T>::max()) {
      return array<uint8_t, 3>{{128, 128, 128}};
    }
    return x < T(0) ?
      array<uint8_t, 3>{{
        uint8_t(0),
        uint8_t(0),
        clamp<uint8_t>(
          T(0),
          T(numeric_limits<uint8_t>::max()),
          numeric_limits<uint8_t>::max() * fabs(x))}} :
      array<uint8_t, 3>{{
        clamp<uint8_t>(
          T(0),
          T(numeric_limits<uint8_t>::max()),
          numeric_limits<uint8_t>::max() * x),
        uint8_t(0),
        uint8_t(0)}};
  };

  auto const width = dist_stats.grid_size[0];
  auto const height = dist_stats.grid_size[1];

  stringstream ss_input;
  ss_input << prefix << "_input_" << typeid(T).name() << ".ppm";
  auto const normalized_input = signedNormalized<T>(
    begin(dist_stats.input_buffer),
    end(dist_stats.input_buffer));
  ppm::writeRgbImage(
    ss_input.str(),
    width,
    height,
    pixelsFromValues(
     begin(normalized_input),
     end(normalized_input),
     pixel_from_value));

  stringstream ss_distance;
  ss_distance << prefix << "_distance_" << typeid(T).name() << ".ppm";
  auto const normalized_distance = signedNormalized<T>(
    begin(dist_stats.distance_buffer),
    end(dist_stats.distance_buffer));
  ppm::writeRgbImage(
    ss_distance.str(),
    width,
    height,
    pixelsFromValues(
      begin(normalized_distance),
      end(normalized_distance),
      pixel_from_value));

  stringstream ss_gt;
  ss_gt << prefix << "_gt_" << typeid(T).name() << ".ppm";
  auto const normalized_gt = signedNormalized<T>(
    begin(dist_stats.distance_ground_truth_buffer),
    end(dist_stats.distance_ground_truth_buffer));
  ppm::writeRgbImage(
    ss_gt.str(),
    width,
    height,
    pixelsFromValues(
      begin(normalized_gt),
      end(normalized_gt),
      pixel_from_value));

  stringstream ss_error;
  ss_error << prefix << "_error_" << typeid(T).name() << ".ppm";
  auto const normalized_error = signedNormalized<T>(
    begin(dist_stats.error_buffer),
    end(dist_stats.error_buffer));
  ppm::writeRgbImage(
    ss_error.str(),
    width,
    height,
    pixelsFromValues(
      begin(normalized_error),
      end(normalized_error),
      pixel_from_value));
}

} // namespace


namespace std {

template<typename T, size_t N>
ostream& operator<<(
  ostream& os,
  thinks::fmm::test::GradientMagnitudeStats<T, N> const& grad_mag_stats)
{
  os << "Gradient magnitude stats <" << typeid(T).name() << ", " << N << ">:" << endl
    << "min abs error: " << grad_mag_stats.min_abs_error << endl
    << "max abs error: " << grad_mag_stats.max_abs_error << endl
    << "avg abs error: " << grad_mag_stats.avg_abs_error << endl
    << "std_dev abs error: " << grad_mag_stats.std_dev_abs_error << endl
    << "duration: " << grad_mag_stats.duration_in_s << " [s]" << endl;
  return os;
}


template<typename T, size_t N>
ostream& operator<<(
  ostream& os,
  thinks::fmm::test::DistanceValueStats<T, N> const& dist_stats)
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


int main(int argc, char* argv[])
{
  using namespace std;
  using namespace thinks::fmm;

  try {
    {
      cout << "Unsigned distance" << endl
           << "-----------------" << endl;

#if 1
      auto const grad_mag_stats2f = test::unsignedGradientMagnitudeStats<float, 2>();
      cout << grad_mag_stats2f << endl;
      writeGradMagImages<float>(grad_mag_stats2f, "unsigned_grad_mag");

      auto const grad_mag_stats2d = test::unsignedGradientMagnitudeStats<double, 2>();
      cout << grad_mag_stats2d << endl;
      writeGradMagImages<double>(grad_mag_stats2d, "unsigned_grad_mag");

      auto const dist_stats2f = test::unsignedDistanceValueStats<float, 2>();
      cout << dist_stats2f << endl;
      writeDistStatImages<float>(dist_stats2f, "unsigned_dist_stat");

      auto const dist_stats2d = test::unsignedDistanceValueStats<double, 2>();
      cout << dist_stats2d << endl;
      writeDistStatImages<double>(dist_stats2d, "unsigned_dist_stat");
#endif

#if 0
      auto const grad_mag_stats3f = test::unsignedGradientMagnitudeStats<float, 3>();
      cout << grad_mag_stats3f << endl;
      auto const grad_mag_stats3d = test::unsignedGradientMagnitudeStats<double, 3>();
      cout << grad_mag_stats3d << endl;

      auto const dist_stats3f = test::unsignedDistanceValueStats<float, 3>();
      cout << dist_stats3f << endl;
      auto const dist_stats3d = test::unsignedDistanceValueStats<double, 3>();
      cout << dist_stats3d << endl;
#endif

#if 0
      auto const grad_mag_stats4f = test::unsignedGradientMagnitudeStats<float, 4>();
      cout << grad_mag_stats4f << endl;
      auto const grad_mag_stats4d = test::unsignedGradientMagnitudeStats<double, 4>();
      cout << grad_mag_stats4d << endl;

      auto const dist_stats4f = test::unsignedDistanceValueStats<float, 4>();
      cout << dist_stats4f << endl;
      auto const dist_stats4d = test::unsignedDistanceValueStats<double, 4>();
      cout << dist_stats4d << endl;
#endif
    }

    {
      cout << "Signed distance" << endl
           << "-----------------" << endl;

#if 1
      auto const grad_mag_stats2f = test::signedGradientMagnitudeStats<float, 2>();
      cout << grad_mag_stats2f << endl;
      writeGradMagImages<float>(grad_mag_stats2f, "signed_grad_mag");

      auto const grad_mag_stats2d = test::signedGradientMagnitudeStats<double, 2>();
      cout << grad_mag_stats2d << endl;
      writeGradMagImages<double>(grad_mag_stats2d, "signed_grad_mag");

      auto const dist_stats2f = test::signedDistanceValueStats<float, 2>();
      cout << dist_stats2f << endl;
      writeDistStatImages<float>(dist_stats2f, "signed_dist_stat");

      auto const dist_stats2d = test::signedDistanceValueStats<double, 2>();
      cout << dist_stats2d << endl;
      writeDistStatImages<double>(dist_stats2d, "signed_dist_stat");
#endif

#if 0
      auto const grad_mag_stats3f = test::signedGradientMagnitudeStats<float, 3>();
      cout << grad_mag_stats3f << endl;
      auto const grad_mag_stats3d = test::signedGradientMagnitudeStats<double, 3>();
      cout << grad_mag_stats3d << endl;

      auto const dist_stats3f = test::signedDistanceValueStats<float, 3>();
      cout << dist_stats3f << endl;
      auto const dist_stats3d = test::signedDistanceValueStats<double, 3>();
      cout << dist_stats3d << endl;
#endif

#if 0
      auto const grad_mag_stats4f = test::signedGradientMagnitudeStats<float, 4>();
      cout << grad_mag_stats4f << endl;
      auto const grad_mag_stats4d = test::signedGradientMagnitudeStats<double, 4>();
      cout << grad_mag_stats4d << endl;

      auto const dist_stats4f = test::signedDistanceValueStats<float, 4>();
      cout << dist_stats4f << endl;
      auto const dist_stats4d = test::signedDistanceValueStats<double, 4>();
      cout << dist_stats4d << endl;
#endif
    }
  }
  catch (exception& ex) {
    cerr << "Exception: " << ex.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
