#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <thinks/testFastMarchingMethod.hpp>


namespace {

template<typename R, typename T> inline
R clamp(T const min, T const max, T const value)
{
  return static_cast<R>(std::min<T>(max, std::max<T>(min, value)));
}


template<typename T, typename R, typename U>
std::vector<R> transformedVector(std::vector<T> const& v, U const unary_op)
{
  using namespace std;

  auto r = vector<R>(v.size());
  transform(begin(v), end(v), begin(r), unary_op);
  return r;
}


template<typename T>
struct Pixel
{
  typedef T ChannelType;
  Pixel() : r(0), g(0), b(0) {}
  Pixel(T const x) : r(x), g(x), b(x) {}
  Pixel(T const _r, T const _g, T const _b) : r(_r), g(_g), b(_b) {}
  T r;
  T g;
  T b;
};

typedef Pixel<std::uint8_t> Pixel8;

void writePpm(std::string const& filename,
              std::size_t const width, std::size_t const height,
              std::vector<Pixel8> const& pixels)
{
  // C-HACK!!
  FILE* fd = fopen(filename.c_str(), "wb");
  (void) fprintf(fd, "P6\n%d %d\n255\n", width, height);
  (void) fwrite((const void*)pixels.data(), sizeof(Pixel8),
                width * height, fd);
  (void) fflush(fd);
  fclose(fd);
}

template<typename T, typename C>
std::vector<Pixel8> pixels(std::vector<T> const& values, C const converter)
{
  using namespace std;

  if (values.empty()) {
    throw runtime_error("empty values vector");
  }

  auto max_value = numeric_limits<T>::lowest();
  auto min_value = numeric_limits<T>::max();
  for (auto iter = begin(values); iter != end(values); ++iter) {
    if (*iter != numeric_limits<T>::max()) {
      max_value = max(max_value, *iter);
      min_value = min(min_value, *iter);
    }
  }

  auto normalized_values = values;
  transform(begin(values), end(values), begin(normalized_values),
    [=](auto const x) {
      if (x == numeric_limits<T>::max()) {
        return x;
      }
      if (x > T(0) && max_value > T(0)) {
        return x / max_value;
      }
      if (x < T(0) && min_value < T(0)) {
        return -(x / min_value);
      }
      return T(0);
    });

  auto pixels = vector<Pixel8>(values.size());
  transform(begin(normalized_values), end(normalized_values), begin(pixels),
            converter);
  return pixels;
}


template<typename T>
void writeGradMagImages(
  thinks::fmm::test::GradientMagnitudeStats<T, 2> const& grad_mag_stats,
  std::string const& prefix)
{
  using namespace std;

  // Negative values in shades of blue, positive values in shades of red.
  auto const pixel_from_value = [](T const x) {
    return x < T(0) ?
      Pixel8(
        Pixel8::ChannelType(0),
        Pixel8::ChannelType(0),
        clamp<Pixel8::ChannelType>(
          T(0),
          T(numeric_limits<Pixel8::ChannelType>::max()),
          numeric_limits<Pixel8::ChannelType>::max() * fabs(x))) :
      Pixel8(
        clamp<Pixel8::ChannelType>(
          T(0),
          T(numeric_limits<Pixel8::ChannelType>::max()),
          numeric_limits<Pixel8::ChannelType>::max() * x),
        Pixel8::ChannelType(0),
        Pixel8::ChannelType(0));
  };

  stringstream ss_input;
  ss_input << prefix << "_input_" << typeid(T).name() << ".ppm";
  writePpm(
    ss_input.str(),
    grad_mag_stats.grid_size[0], grad_mag_stats.grid_size[1],
    pixels(
      grad_mag_stats.input_buffer,
      [=](T const d) {
        if (d == numeric_limits<T>::max()) {
          return Pixel8(
            Pixel8::ChannelType(128),
            Pixel8::ChannelType(128),
            Pixel8::ChannelType(128));
        }
        return pixel_from_value(d);
    }));

  stringstream ss_distance;
  ss_distance << prefix << "_distance_" << typeid(T).name() << ".ppm";
  writePpm(
    ss_distance.str(),
    grad_mag_stats.grid_size[0], grad_mag_stats.grid_size[1],
    pixels(grad_mag_stats.distance_buffer, pixel_from_value));

  stringstream ss_grad_mag;
  ss_grad_mag << prefix << "_" << typeid(T).name() << ".ppm";
  writePpm(
    ss_grad_mag.str(),
    grad_mag_stats.grid_size[0], grad_mag_stats.grid_size[1],
    pixels(
      transformedVector<array<T, 2>, T>(
        grad_mag_stats.grad_buffer,
          [](auto const v) { return sqrt(v[0] * v[0] + v[1] * v[1]); } ),
      [](T const& v) {
        return Pixel8(
          clamp<Pixel8::ChannelType>(
            T(0),
            T(numeric_limits<Pixel8::ChannelType>::max()),
            numeric_limits<Pixel8::ChannelType>::max() * v));
      }));

  stringstream ss_error;
  ss_error << prefix << "_error_" << typeid(T).name() << ".ppm";
  writePpm(
    ss_error.str(),
    grad_mag_stats.grid_size[0], grad_mag_stats.grid_size[1],
    pixels(grad_mag_stats.error_buffer, pixel_from_value));
}

template<typename T>
void writeDistStatImages(
  thinks::fmm::test::DistanceValueStats<T, 2> const& dist_stats,
  std::string const& prefix)
{
  using namespace std;

  auto const pixel_from_distance = [](T const d) {
    return d < T(0) ?
      Pixel8(
        Pixel8::ChannelType(0),
        Pixel8::ChannelType(0),
        clamp<Pixel8::ChannelType>(
          T(0),
          T(numeric_limits<Pixel8::ChannelType>::max()),
          numeric_limits<Pixel8::ChannelType>::max() * fabs(d))) :
      Pixel8(
        clamp<Pixel8::ChannelType>(
          T(0),
          T(numeric_limits<Pixel8::ChannelType>::max()),
          numeric_limits<Pixel8::ChannelType>::max() * d),
        Pixel8::ChannelType(0),
        Pixel8::ChannelType(0));
  };

  stringstream ss_input;
  ss_input << prefix << "_input_" << typeid(T).name() << ".ppm";
  writePpm(
    ss_input.str(),
    dist_stats.grid_size[0], dist_stats.grid_size[1],
    pixels(
      dist_stats.input_buffer,
      [=](T const d) {
        if (d == numeric_limits<T>::max()) {
          return Pixel8(
            Pixel8::ChannelType(128),
            Pixel8::ChannelType(128),
            Pixel8::ChannelType(128));
        }
        return pixel_from_distance(d);
      }));

  stringstream ss_distance;
  ss_distance << prefix << "_distance_" << typeid(T).name() << ".ppm";
  writePpm(
    ss_distance.str(),
    dist_stats.grid_size[0], dist_stats.grid_size[1],
    pixels(dist_stats.distance_buffer, pixel_from_distance));

  stringstream ss_gt;
  ss_gt << prefix << "_gt_" << typeid(T).name() << ".ppm";
  writePpm(
    ss_gt.str(),
    dist_stats.grid_size[0], dist_stats.grid_size[1],
    pixels(dist_stats.distance_ground_truth_buffer, pixel_from_distance));

  stringstream ss_error;
  ss_error << prefix << "_error_" << typeid(T).name() << ".ppm";
  writePpm(
    ss_error.str(),
    dist_stats.grid_size[0], dist_stats.grid_size[1],
    pixels(dist_stats.error_buffer, pixel_from_distance));
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
  using namespace thinks::fmm::test;

  try {
#if 1
    {
      cout << "Unsigned distance" << endl
           << "-----------------" << endl;

      auto const grad_mag_stats2f = UnsignedGradientMagnitudeStats<float, 2>();
      cout << grad_mag_stats2f << endl;
      writeGradMagImages<float>(grad_mag_stats2f, "unsigned_grad_mag");

      auto const grad_mag_stats2d = UnsignedGradientMagnitudeStats<double, 2>();
      cout << grad_mag_stats2d << endl;
      writeGradMagImages<double>(grad_mag_stats2d, "unsigned_grad_mag");

#if 0
      auto const grad_mag_stats3f = UnsignedGradientMagnitudeStats<float, 3>();
      cout << grad_mag_stats3f << endl;
      auto const grad_mag_stats3d = UnsignedGradientMagnitudeStats<double, 3>();
      cout << grad_mag_stats3d << endl;
#endif

#if 0
      Gradient magnitude stats <double, 3>:
      min: 0.0204039
      max: 1.02286
      avg: 0.9975
      std_dev: 0.0364526
      duration: 5.42531 [s]


      auto const grad_mag_stats4f = UnsignedGradientMagnitudeStats<float, 4>();
      cout << grad_mag_stats4f << endl;
      auto const grad_mag_stats4d = UnsignedGradientMagnitudeStats<double, 4>();
      cout << grad_mag_stats4d << endl;
#endif

      auto const dist_stats2f = UnsignedDistanceValueStats<float, 2>();
      cout << dist_stats2f << endl;
      writeDistStatImages<float>(dist_stats2f, "unsigned_dist_stat");

      auto const dist_stats2d = UnsignedDistanceValueStats<double, 2>();
      cout << dist_stats2d << endl;
      writeDistStatImages<double>(dist_stats2d, "unsigned_dist_stat");

#if 0
      auto const dist_stats3f = UnsignedDistanceValueStats<float, 3>();
      cout << dist_stats3f << endl;
      auto const dist_stats3d = UnsignedDistanceValueStats<double, 3>();
      cout << dist_stats3d << endl;
#endif

#if 0
      auto const dist_stats4f = UnsignedDistanceValueStats<float, 4>();
      cout << dist_stats4f << endl;
      auto const dist_stats4d = UnsignedDistanceValueStats<double, 4>();
      cout << dist_stats4d << endl;
#endif
    }
#endif
    {
      cout << "Signed distance" << endl
           << "-----------------" << endl;

      auto const grad_mag_stats2f = SignedGradientMagnitudeStats<float, 2>();
      cout << grad_mag_stats2f << endl;
      writeGradMagImages<float>(grad_mag_stats2f, "signed_grad_mag");

      auto const grad_mag_stats2d = SignedGradientMagnitudeStats<double, 2>();
      cout << grad_mag_stats2d << endl;
      writeGradMagImages<double>(grad_mag_stats2d, "signed_grad_mag");

#if 0
      auto const grad_mag_stats3f = SignedGradientMagnitudeStats<float, 3>();
      cout << grad_mag_stats3f << endl;
      auto const grad_mag_stats3d = SignedGradientMagnitudeStats<double, 3>();
      cout << grad_mag_stats3d << endl;
#endif
      //thinks::fmm::test::testUnsignedGradientMagnitude<float, 4>();
      //thinks::fmm::test::testUnsignedGradientMagnitude<double, 4>();

      auto const dist_stats2f = SignedDistanceValueStats<float, 2>();
      cout << dist_stats2f << endl;
      writeDistStatImages<float>(dist_stats2f, "signed_dist_stat");

      auto const dist_stats2d = SignedDistanceValueStats<double, 2>();
      cout << dist_stats2d << endl;
      writeDistStatImages<double>(dist_stats2d, "signed_dist_stat");

#if 0
      auto const dist_stats3f = SignedDistanceValueStats<float, 3>();
      cout << dist_stats3f << endl;
      auto const dist_stats3d = SignedDistanceValueStats<double, 3>();
      cout << dist_stats3d << endl;
#endif

      //thinks::fmm::test::testUnsignedDistanceValues<float, 4>();
      //thinks::fmm::test::testUnsignedDistanceValues<double, 4>();
    }
  }
  catch (exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }

  return 0;
}


