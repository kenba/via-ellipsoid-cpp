//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2019-2024 Ken Barker
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
/// @brief Contains accuracy tests for ellipsoid::calculate_azimuth_aux_length.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid/geodesic_functions.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/test/unit_test.hpp>
#include <filesystem>
#include <fstream>

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
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_geodesic_examples_suite)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_examples) {
  // Test reads file from directory, fails if directory environment variable is
  // not present.
  const auto env_var(std::getenv("GEODTEST_DIR"));
  BOOST_REQUIRE(env_var != nullptr);

  const auto filename("GeodTest.dat");
  std::filesystem::path file_path(env_var);
  file_path /= filename;

  std::ifstream data_file(file_path.c_str());
  BOOST_REQUIRE(data_file.is_open());

  const auto WGS84_ELLIPSOID{ellipsoid::Ellipsoid<double>::wgs84()};

  auto line_number(0u);
  while (!data_file.eof()) {
    // read next line from the file
    std::string line;
    std::getline(data_file, line);
    if (!line.empty()) {
      const auto params{tokenise(line)};
      const double lat1{std::stod(params[LAT_1])};
      const double lon1{std::stod(params[LON_1])};
      const double azimuth_1{std::stod(params[AZI_1])};

      const double lat2{std::stod(params[LAT_2])};
      const double lon2{std::stod(params[LON_2])};

      const double distance_m{std::stod(params[D_METRES])};
      const double distance_deg{std::stod(params[D_DEGREES])};

      const Degrees<double> lat1d{lat1};
      const Degrees<double> lon1d{lon1};
      const Degrees<double> lat2d{lat2};
      const Degrees<double> lon2d{lon2};
      const auto [azimuth, aux_length,
                  iterations]{ellipsoid::calculate_azimuth_aux_length(
          LatLong(lat1d, lon1d), LatLong(lat2d, lon2d), WGS84_ELLIPSOID)};

      // Compare azimuth
      const double delta_azimuth{
          std::abs(azimuth_1 - azimuth.to_degrees().v())};
      const double azimuth_tolerance{(line_number <= 400000) ? 5.331e-5
                                                             : 0.077};
      BOOST_CHECK_SMALL(delta_azimuth, 100 * azimuth_tolerance);

      // Compare aux sphere great circle length
      const double delta_aux_length{
          std::abs(trig::deg2rad(distance_deg) - aux_length.v())};
      const double aux_length_tolerance{3.0e-10};
      BOOST_CHECK_SMALL(delta_aux_length, 100 * aux_length_tolerance);

      // Compare geodesic length
      const auto beta1{
          WGS84_ELLIPSOID.calculate_parametric_latitude(Angle(lat1d))};
      const auto result_m{ellipsoid::convert_radians_to_metres(
          beta1, azimuth, aux_length, WGS84_ELLIPSOID)};
      const double delta_length_m{std::abs(distance_m - result_m.v())};

      // if a short geodesic, compare delta length
      if (line_number >= 150000 && line_number < 200000) {
        BOOST_CHECK_SMALL(delta_length_m, 100 * 9.0e-5);
      } else {
        BOOST_CHECK_CLOSE(distance_m, result_m.v(), 100 * 2.5e-9);
      }
    }

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

  std::cout << "lines: " << line_number << std::endl;
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
