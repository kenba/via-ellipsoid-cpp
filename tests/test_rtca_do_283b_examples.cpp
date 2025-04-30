//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024-2025 Ken Barker
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
/// @file test_rtca_do_283b_examples.cpp
/// @brief Contains geodesic accuracy tests using RTCA DO-283B Appendix E
/// examples.
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
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_rtca_do_283b_geodesic_examples_suite)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_rtca_do_283b_geodesic_examples) {
  const auto env_var(std::getenv("ViaEllipsoid_DIR"));
  BOOST_REQUIRE(env_var != nullptr);

  const auto filename{"rtca_do_283b_geodesic_examples.csv"};
  std::filesystem::path file_path(env_var);
  file_path /= "data";
  file_path /= filename;

  std::ifstream data_file(file_path.c_str());
  BOOST_REQUIRE(data_file.is_open());

  const auto WGS84_ELLIPSOID{ellipsoid::Ellipsoid<double>::wgs84()};

  std::string line;
  // ignore first line in the file
  std::getline(data_file, line);
  while (!data_file.eof()) {
    // read next line from the file
    std::getline(data_file, line);
    if (!line.empty()) {
      const auto params{tokenise(line)};
      const double lat1{std::stod(params[0])};
      const double lon1{std::stod(params[1])};
      const double lat2{std::stod(params[2])};
      const double lon2{std::stod(params[3])};
      const double azimuth_1{std::stod(params[4])};
      // const double azimuth_2{std::stod(params[5])};
      const double d_metres{std::stod(params[6])};
      // const double d_error{std::stod(params[7])};

      const Degrees<double> lat1d{lat1};
      const Degrees<double> lon1d{lon1};
      const Degrees<double> lat2d{lat2};
      const Degrees<double> lon2d{lon2};
      const auto [azimuth, aux_length, _end_azimuth, iterations]{
          ellipsoid::calculate_azimuths_arc_length(LatLong(lat1d, lon1d),
                                                   LatLong(lat2d, lon2d))};

      // Compare azimuth
      const double delta_azimuth{
          std::abs(azimuth_1 - azimuth.to_degrees().v())};
      const double azimuth_tolerance{1e-9};
      BOOST_CHECK_SMALL(delta_azimuth, 100 * azimuth_tolerance);

      // Compare geodesic length
      const auto beta1{
          WGS84_ELLIPSOID.calculate_parametric_latitude(Angle(lat1d))};
      const auto result_m{
          ellipsoid::convert_radians_to_metres(beta1, azimuth, aux_length)};
      const double delta_length_m{std::abs(d_metres - result_m.v())};

      const double distance_tolerance{1e-5};
      BOOST_CHECK_SMALL(delta_length_m, 100 * distance_tolerance);
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
