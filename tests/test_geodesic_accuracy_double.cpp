//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2019-2025 Ken Barker
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
/// @file test_geodesic_accuracy_double.cpp
/// @brief Contains accuracy tests for ellipsoid::calculate_azimuths_arc_length.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid/geodesic_functions.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/test/unit_test.hpp>
#include <filesystem>
#include <fstream>
#include <via/angle/trig.hpp>
#if defined(OUTPUT_VIA_ELLIPSOID_VALUES) || defined(OUTPUT_GEOGRAPHICLIB_VALUES)
#include <iomanip>
#endif
#ifdef TEST_VINCENTY
#include "via/ellipsoid/vincenty_functions.hpp"
#endif
#ifdef TEST_GEOGRAPHICLIB
#include <GeographicLib/Geodesic.hpp>
#endif

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

using PositionData = std::array<double, 8>;

/// Read PositionData from the first 10000 lines of data_file.
std::vector<PositionData> read_position_data(std::ifstream &data_file) {
  std::vector<PositionData> data{};

  const auto line_count(500'000u);
  data.reserve(line_count);

  // Read data_file into data
  while (!data_file.eof()) {
    // read next line from the file
    std::string line;
    std::getline(data_file, line);
    if (!line.empty()) {
      const auto params{tokenise(line)};
      // Start position and azimuth
      const double lat1{std::stod(params[LAT_1])};
      const double lon1{std::stod(params[LON_1])};
      const double azimuth_1{std::stod(params[AZI_1])};

      // End position and azimuth
      const double lat2{std::stod(params[LAT_2])};
      const double lon2{std::stod(params[LON_2])};
      const double azimuth_2{std::stod(params[AZI_2])};

      const double distance_m{std::stod(params[D_METRES])};
      const double distance_deg{std::stod(params[D_DEGREES])};

      data.emplace_back(PositionData{lat1, lon1, azimuth_1, lat2, lon2,
                                     azimuth_2, distance_m, distance_deg});
    }
  }

  return data;
}
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_geodesic_examples_suite)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_examples) {
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

  auto line_number(0u);

  for (const auto &position : data) {
    const Degrees<double> lat1d{position[LAT_1]};
    const Degrees<double> lon1d{position[LON_1]};
    const Degrees<double> lat2d{position[LAT_2]};
    const Degrees<double> lon2d{position[LON_2]};

    const auto [azimuth, aux_length, end_azimuth, iterations]{
        ellipsoid::calculate_azimuths_arc_length(LatLong(lat1d, lon1d),
                                                 LatLong(lat2d, lon2d))};

    // Compare azimuth
    const double azimuth_1{position[AZI_1]};
    const double delta_azimuth{std::abs(azimuth_1 - azimuth.to_degrees().v())};
#ifdef _MSC_VER
    // relax tolerance for Visual Studio
    const double azimuth_tolerance{(line_number <= 400000) ? 5.33e-7 : 3.1e-4};
#else
    const double azimuth_tolerance{(line_number <= 400000) ? 5.0e-7 : 3.1e-4};
#endif
    BOOST_CHECK_SMALL(delta_azimuth, 100 * azimuth_tolerance);

    // Compare aux sphere great circle length
    const double distance_deg{position[D_DEGREES]};
    const double delta_aux_length{
        std::abs(trig::deg2rad(distance_deg) - aux_length.v())};
    const double aux_length_tolerance{3.0e-10};
    BOOST_CHECK_SMALL(delta_aux_length, 100 * aux_length_tolerance);

    // Compare geodesic length
    const double distance_m{position[D_METRES]};
    const auto beta1{
        WGS84_ELLIPSOID.calculate_parametric_latitude(Angle(lat1d))};
    const auto result_m{
        ellipsoid::convert_radians_to_metres(beta1, azimuth, aux_length)};
    const double delta_length_m{std::abs(distance_m - result_m.v())};

    // if a short geodesic, compare delta length
    if (line_number >= 150000 && line_number < 200000) {
      BOOST_CHECK_SMALL(delta_length_m, 100 * 3.4e-11);
    } else {
#ifdef _MSC_VER
      // relax tolerance for Visual Studio
      BOOST_CHECK_CLOSE(distance_m, result_m.v(), 100 * 5.9e-13);
#else
      BOOST_CHECK_CLOSE(distance_m, result_m.v(), 100 * 5.0e-13);
#endif
    }

    // Compare end azimuth
    const double azimuth_2{position[AZI_2]};
    const double delta_end_azimuth{
        std::abs(azimuth_2 - end_azimuth.to_degrees().v())};
    BOOST_CHECK_SMALL(delta_end_azimuth, 100 * azimuth_tolerance);

#ifdef OUTPUT_VIA_ELLIPSOID_VALUES
    std::cout << std::setprecision(18) << lat1d << " 0 "
              << azimuth.to_degrees().v() << " " << lat2d << " " << lon2d << " "
              << end_azimuth.to_degrees().v() << " " << result_m.v()
              << std::endl;
#endif

    //  random_df = tests_df[:100000]
    //  antipodal_df = tests_df[100000:150000]
    //  short_df = tests_df[150000:200000]
    //  one_pole_df = tests_df[200000:250000]
    //  two_poles_df = tests_df[250000:300000]
    //  near_meridional_df = tests_df[300000:350000]
    //  near_equatorial_df = tests_df[350000:400000]
    //  between_vertices_df = tests_df[400000:450000]
    //  end_by_vertices_df = tests_df[450000:500000]
    ++line_number;
  }

#ifndef OUTPUT_VIA_ELLIPSOID_VALUES
#ifndef OUTPUT_GEOGRAPHICLIB_VALUES
  std::cout << "lines: " << line_number << std::endl;
#endif
#endif
}
//////////////////////////////////////////////////////////////////////////////

#ifdef TEST_VINCENTY
//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_examples_with_boost_vincenty) {
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

  auto line_number(0u);

  for (const auto &position : data) {
    const Degrees<double> lat1d{position[LAT_1]};
    const Degrees<double> lon1d{position[LON_1]};
    const Degrees<double> lat2d{position[LAT_2]};
    const Degrees<double> lon2d{position[LON_2]};

    const auto [azimuth_rad, s12, _rev_azimuth]{
        ellipsoid::vincenty::inverse_azimuths_and_distance(
            LatLong<double>(lat1d, lon1d), LatLong<double>(lat2d, lon2d))};

    // Compare azimuth
    const double azimuth_1{position[AZI_1]};
    const double delta_azimuth{
        std::abs(azimuth_1 - trig::rad2deg(azimuth_rad.v()))};
    const double azimuth_tolerance{5.1};
    BOOST_CHECK_SMALL(delta_azimuth, 100 * azimuth_tolerance);

    // Compare geodesic length
    const double distance_m{position[D_METRES]};
    const double delta_length_m{distance_m - s12.v()};

    // if a short geodesic, compare delta length
    if (line_number >= 150000 && line_number < 200000) {
      BOOST_CHECK_SMALL(delta_length_m, 100 * 1.0e-7);
    } else {
      BOOST_CHECK_CLOSE(distance_m, s12.v(), 100 * 4.0e-3);
    }

    ++line_number;
  }
}
//////////////////////////////////////////////////////////////////////////////
#endif

#ifdef TEST_GEOGRAPHICLIB
//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_examples_with_geographiclib) {
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

  auto line_number(0u);

  for (const auto &position : data) {
    const double lat1d{position[LAT_1]};
    const double lon1d{position[LON_1]};
    const double azimuth_1{position[AZI_1]};
    const double lat2d{position[LAT_2]};
    const double lon2d{position[LON_2]};

    // Solve the inverse problem on the WGS-84 geoid using GeographicLib
    double s12, azi1, azi2;
    geoid.Inverse(lat1d, lon1d, lat2d, lon2d, s12, azi1, azi2);

    // Compare geodesic length
    const double distance_m{position[D_METRES]};
    const double delta_length_m{distance_m - s12};

    // if a short geodesic, compare delta length
    if (line_number >= 150000 && line_number < 200000) {
      BOOST_CHECK_SMALL(delta_length_m, 100 * 3.0e-11);
    } else {
      BOOST_CHECK_CLOSE(distance_m, s12, 100 * 5.3e-13);
    }

    // Compare start azimuth
    const double delta_azimuth1{azimuth_1 - azi1};
    const double azimuth_tolerance{(line_number <= 400000) ? 5.0e-7 : 2.0e-4};
    BOOST_CHECK_SMALL(delta_azimuth1, 100 * azimuth_tolerance);

    // Compare end azimuth
    const double azimuth_2{position[AZI_2]};
    const double delta_azimuth2{azimuth_2 - azi2};
    BOOST_CHECK_SMALL(delta_azimuth2, 100 * azimuth_tolerance);

#ifdef OUTPUT_GEOGRAPHICLIB_VALUES
    std::cout << std::setprecision(18) << lat1d << " 0 " << azi1 << " " << lat2d
              << " " << lon2d << " " << azi2 << " " << s12 << std::endl;
#endif

#ifdef OUTPUT_QUADRANT_MISMATCHES
    // Output lines where GeographicLib and GeodTest.dat start azimuth
    // quadrants are different
    if ((azimuth_1 > 90 && azi1 < 90) || (azimuth_1 < 90 && azi1 > 90)) {
      std::cout << "Line: " << line_number << std::setprecision(15)
                << " delta length: " << delta_length_m
                << " delta azimuth: " << delta_azimuth1 << std::fixed
                << " ref azimuth: " << azimuth_1 << " azimuth: " << azi1
                << " ref end azimuth: " << azimuth_2 << " end azimuth: " << azi2
                << " lat1: " << lat1d << " lat2: " << lat2d << " lon: " << lon2d
                << std::endl;
    }
#endif

    ++line_number;
  }

#ifndef OUTPUT_GEOGRAPHICLIB_VALUES
  std::cout << "GeographicLib lines: " << line_number << std::endl;
#endif
}
//////////////////////////////////////////////////////////////////////////////
#endif

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
