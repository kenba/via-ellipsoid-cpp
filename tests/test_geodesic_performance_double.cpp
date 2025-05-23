//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2025 Ken Barker
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//////////////////////////////////////////////////////////////////////////////
/// @file test_geodesic_performance_double.cpp
/// @brief Contains performance tests for geodesics.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid.hpp"
#include <array>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/test/unit_test.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <vector>
#ifdef TEST_VINCENTY
#include "via/ellipsoid/vincenty_functions.hpp"
#endif
#ifdef TEST_GEOGRAPHICLIB
#include <GeographicLib/Intersect.hpp>
#endif

using namespace std::chrono;
using namespace via;

namespace {
// A function to tokenise a line from the data file
std::vector<std::string> tokenise(std::string const &s) {
  std::vector<std::string> tokens;
  boost::split(tokens, s, boost::is_any_of(", "), boost::token_compress_on);
  return tokens;
}

// The columns of the data file.
enum Columns {
  LAT_1,
  LON_1,
  AZI_1,
  LAT_2,
  LON_2,
  AZI_2,
  D_METRES,
  D_DEGREES,
  M12,
  AREA
};

using PositionData = std::array<double, 5>;

/// Read PositionData from the first 10000 lines of data_file.
std::vector<PositionData> read_position_data(std::ifstream &data_file) {
  std::vector<PositionData> data{};

  const auto random_line_count(100'000u);
  data.reserve(random_line_count);

  // Read data_file into data
  auto line_number(0u);
  while (!data_file.eof() && line_number < random_line_count) {
    // read next line from the file
    std::string line;
    std::getline(data_file, line);
    if (!line.empty()) {
      const auto params{tokenise(line)};
      // Start position and azimuth
      const double lat1{std::stod(params[LAT_1])};
      const double lon1{std::stod(params[LON_1])};
      const double azimuth_1{std::stod(params[AZI_1])};

      // End position
      const double lat2{std::stod(params[LAT_2])};
      const double lon2{std::stod(params[LON_2])};
      data.emplace_back(PositionData{lat1, lon1, azimuth_1, lat2, lon2});
    }

    ++line_number;
  }

  return data;
}
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_geodesic_examples_suite)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_intersection_performance) {
  // Test reads file from directory, fails if directory environment variable is
  // not present.
  const auto env_var(std::getenv("GEODTEST_DIR"));
  BOOST_REQUIRE(env_var != nullptr);

  const auto filename("GeodTest.dat");
  std::filesystem::path file_path(env_var);
  file_path /= filename;

  std::ifstream data_file(file_path.c_str());
  BOOST_REQUIRE(data_file.is_open());

  // The position data
  std::vector<PositionData> data{read_position_data(data_file)};

  // The geodesic data
  std::vector<ellipsoid::GeodesicSegment<double>> geodesic_data{};
  geodesic_data.reserve(data.size());

  // Create random Geodesics
  const auto t0{high_resolution_clock::now()};
  for (const auto &position : data) {
    const Degrees<double> lat1d{position[LAT_1]};
    const Degrees<double> lon1d{position[LON_1]};
    const Degrees<double> lat2d{position[LAT_2]};
    const Degrees<double> lon2d{position[LON_2]};
    geodesic_data.emplace_back(ellipsoid::GeodesicSegment<double>(
        LatLong(lat1d, lon1d), LatLong(lat2d, lon2d)));
  }
  const auto t1{high_resolution_clock::now()};

  // calculate time to create random Geodesics
  const duration<double, std::milli> create_geodesics_ms{t1 - t0};
  const auto average_create_geodesics_us{create_geodesics_ms.count() * 1e3 /
                                         geodesic_data.size()};
  std::cout << "Geodesic_data size: " << geodesic_data.size()
            << " time: " << create_geodesics_ms.count() << " ms\n"
            << "Average time per GeodesicSegment: "
            << average_create_geodesics_us << " us" << std::endl;

  // Create a reference GeodesicSegment
  // The longest geodesic from the DO-238B set translated to start at 90W
  const LatLong a(Degrees(1.0), Degrees(-90.0));
  const LatLong b(Degrees(-0.998286322222), Degrees(89.296674991667));
  const ellipsoid::GeodesicSegment<double> reference(a, b);

  // 1mm precision in Radians on the auxiliary sphere
  const units::si::Metres<double> precision_1mm{1e-3};
  // Create intersections with random Geodesics at 1mm precision
  const auto t2{high_resolution_clock::now()};
  for (const auto &geodesic : geodesic_data) {
    const auto result{ellipsoid::calculate_intersection_point(
        reference, geodesic, precision_1mm)};
  }
  const auto t3{high_resolution_clock::now()};

  // calculate time to calculate geodesic intersection points
  duration<double, std::milli> geodesic_intersections_ms{t3 - t2};
  auto average_geodesic_intersections_us{geodesic_intersections_ms.count() *
                                         1e3 / geodesic_data.size()};
  std::cout << "Intersections time, 1mm precision: "
            << geodesic_intersections_ms.count() << " ms\n"
            << "Average time per intersection: "
            << average_geodesic_intersections_us << " us" << std::endl;

  // 1m precision in Radians on the auxiliary sphere
  const units::si::Metres<double> precision_1m{1.0};
  // Create intersections with random Geodesics at 1mm precision
  const auto t4{high_resolution_clock::now()};
  for (const auto &geodesic : geodesic_data) {
    const auto result{ellipsoid::calculate_intersection_point(
        reference, geodesic, precision_1m)};
  }
  const auto t5{high_resolution_clock::now()};

  // calculate time to calculate geodesic intersection points
  geodesic_intersections_ms = t5 - t4;
  average_geodesic_intersections_us =
      geodesic_intersections_ms.count() * 1e3 / geodesic_data.size();
  std::cout << "Intersections time, 1m precision: "
            << geodesic_intersections_ms.count() << " ms\n"
            << "Average time per intersection: "
            << average_geodesic_intersections_us << " us\n"
            << std::endl;
}
//////////////////////////////////////////////////////////////////////////////

#ifdef TEST_VINCENTY
//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_vincenty_inverse_performance) {
  const auto WGS84_ELLIPSOID{ellipsoid::Ellipsoid<double>::wgs84()};

  // Test reads file from directory, fails if directory environment variable is
  // not present.
  const auto env_var(std::getenv("GEODTEST_DIR"));
  BOOST_REQUIRE(env_var != nullptr);

  const auto filename("GeodTest.dat");
  std::filesystem::path file_path(env_var);
  file_path /= filename;

  std::ifstream data_file(file_path.c_str());
  BOOST_REQUIRE(data_file.is_open());

  // The position data
  std::vector<PositionData> data{read_position_data(data_file)};

  // The geodesic data
  std::vector<double> geodesic_data{};
  geodesic_data.reserve(data.size());

  // Create random Geodesics
  const auto t0{high_resolution_clock::now()};
  for (const auto &position : data) {
    const Degrees<double> lat1d{position[LAT_1]};
    const Degrees<double> lon1d{position[LON_1]};
    const Degrees<double> lat2d{position[LAT_2]};
    const Degrees<double> lon2d{position[LON_2]};

    const auto distance{ellipsoid::vincenty::inverse_distance(
        LatLong(lat1d, lon1d), LatLong(lat2d, lon2d))};
    geodesic_data.emplace_back(distance.v());
  }
  const auto t1{high_resolution_clock::now()};

  // calculate time to create random Geodesics
  const duration<double, std::milli> create_geodesics_ms{t1 - t0};
  const auto average_create_geodesics_us{create_geodesics_ms.count() * 1e3 /
                                         geodesic_data.size()};
  std::cout << "vincenty geodesic_data size: " << geodesic_data.size()
            << " time: " << create_geodesics_ms.count() << " ms\n"
            << "vincenty average time per GeodesicLine: "
            << average_create_geodesics_us << " us" << std::endl;
}
//////////////////////////////////////////////////////////////////////////////
#endif

#ifdef TEST_GEOGRAPHICLIB
//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geographiclib_intersection_performance) {
  // The GeographicLib WGS-84 ellipsoid
  const auto geoid(GeographicLib::Geodesic::WGS84());
  // Test reads file from directory, fails if directory environment variable is
  // not present.
  const auto env_var(std::getenv("GEODTEST_DIR"));
  BOOST_REQUIRE(env_var != nullptr);

  const auto filename("GeodTest.dat");
  std::filesystem::path file_path(env_var);
  file_path /= filename;

  std::ifstream data_file(file_path.c_str());
  BOOST_REQUIRE(data_file.is_open());

  // The position data
  std::vector<PositionData> data{read_position_data(data_file)};

  // The geodesic data
  std::vector<GeographicLib::GeodesicLine> geodesic_data{};
  geodesic_data.reserve(data.size());

  // Create random Geodesics
  const auto t0{high_resolution_clock::now()};
  for (const auto &position : data) {
    const double lat1d{position[LAT_1]};
    const double lon1d{position[LON_1]};
    const double lat2d{position[LAT_2]};
    const double lon2d{position[LON_2]};
    geodesic_data.emplace_back(geoid.InverseLine(
        lat1d, lon1d, lat2d, lon2d, GeographicLib::Intersect::LineCaps));
  }
  const auto t1{high_resolution_clock::now()};

  // calculate time to create random Geodesics
  const duration<double, std::milli> create_geodesics_ms{t1 - t0};
  const auto average_create_geodesics_us{create_geodesics_ms.count() * 1e3 /
                                         geodesic_data.size()};
  std::cout << "GeographicLib geodesic_data size: " << geodesic_data.size()
            << " time: " << create_geodesics_ms.count() << " ms\n"
            << "GeographicLib average time per GeodesicLine: "
            << average_create_geodesics_us << " us" << std::endl;

  // Create a reference GeodesicLine
  // The longest geodesic from the DO-238B set translated to start at 90W
  const double lat1d{1.0};
  const double lon1d{-90.0};
  const double lat2d{-0.998286322222};
  const double lon2d{89.296674991667};
  const GeographicLib::GeodesicLine reference(geoid.InverseLine(
      lat1d, lon1d, lat2d, lon2d, GeographicLib::Intersect::LineCaps));

  GeographicLib::Intersect intersect(geoid);

  // Create intersections with random GeodesicLines
  const auto t2{high_resolution_clock::now()};
  for (const auto &geodesic : geodesic_data) {
    const auto _point{intersect.Closest(reference, geodesic)};
  }
  const auto t3{high_resolution_clock::now()};

  // calculate time to calculate geodesic intersection points
  const duration<double, std::milli> geodesic_intersections_ms{t3 - t2};
  const auto average_geodesic_intersections_us{
      geodesic_intersections_ms.count() * 1e3 / geodesic_data.size()};
  std::cout << "GeographicLib intersections time: "
            << geodesic_intersections_ms.count() << " ms\n"
            << "GeographicLib average time per intersection: "
            << average_geodesic_intersections_us << " us" << std::endl;
}
//////////////////////////////////////////////////////////////////////////////
#endif

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
