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
/// @file test_ellipsoid_double.cpp
/// @brief Contains unit tests for the via::ellipsoid Geodesic class.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid.hpp"
#include <boost/test/unit_test.hpp>

using namespace via::ellipsoid;
using namespace via;

namespace {
const auto EPSILON(std::numeric_limits<double>::epsilon());
const auto CALCULATION_TOLERANCE(100 * EPSILON);
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_ellipsoid_double)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_atd_and_xtd_karney) {
  // Karney's example
  // Istanbul, Washington and Reyjavik
  const LatLong istanbul(Degrees(42.0), Degrees(29.0));
  const LatLong washington(Degrees(39.0), Degrees(-77.0));
  const Geodesic<double> g1(istanbul, washington);
  BOOST_CHECK(g1.is_valid());

  const LatLong reyjavik(Degrees(64.0), Degrees(-22.0));

  // Calculate geodesic along track and across track distances to 1mm precision.
  const auto [atd, xtd, iterations]{
      g1.calculate_atd_and_xtd(reyjavik, units::si::Metres(1e-3))};
  BOOST_CHECK_CLOSE(3928788.572, atd.v(), 100 * 1e-3);
  BOOST_CHECK_CLOSE(-1010585.9988368, xtd.v(), 100 * 1e-3);

  // Karney's latitude and longitude from Final result at:
  // https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#8a93
  const auto position{g1.lat_long(atd)};
  BOOST_CHECK_CLOSE(54.92853149711691, position.lat().v(),
                    2 * CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(-21.93729106604878, position.lon().v(),
                    200 * CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_closest_intersection_lengths_baselga_1) {
  // First example from Baselga paper
  const LatLong a1(Degrees(52.0), Degrees(5.0));
  const LatLong b1(Degrees(51.4), Degrees(6.0));
  const LatLong a2(Degrees(51.5), Degrees(4.5));
  const LatLong b2(Degrees(52.0), Degrees(5.5));

  const Geodesic<double> g1(a1, b1);
  const Geodesic<double> g2(a2, b2);

  const auto result1{
      calculate_intersection_point(g1, g2, units::si::Metres(1e-3))};
  BOOST_CHECK(result1.has_value());

  BOOST_CHECK_CLOSE(51.86566538889, result1.value().lat().v(), 100 * 1e-3);
  BOOST_CHECK_CLOSE(5.22745711111, result1.value().lon().v(), 100 * 1e-3);

  // Swap geodesics
  const auto result2{
      calculate_intersection_point(g2, g1, units::si::Metres(1e-3))};
  BOOST_CHECK(result2.has_value());

  BOOST_CHECK_CLOSE(51.86566538889, result2.value().lat().v(), 100 * 1e-3);
  BOOST_CHECK_CLOSE(5.22745711111, result2.value().lon().v(), 100 * 1e-3);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_closest_intersection_point_karney) {
  // Second example from Baselga paper
  // Karney's example
  // Istanbul, Washington, Reyjavik and Accra
  const LatLong istanbul(Degrees(42.0), Degrees(29.0));
  const LatLong washington(Degrees(39.0), Degrees(-77.0));
  const LatLong reyjavik(Degrees(64.0), Degrees(-22.0));
  const LatLong accra(Degrees(6.0), Degrees(0.0));

  const Geodesic<double> g1(istanbul, washington);
  const Geodesic<double> g2(reyjavik, accra);

  const auto result1{
      calculate_intersection_point(g1, g2, units::si::Metres(1e-3))};
  BOOST_CHECK(result1.has_value());

  BOOST_CHECK_CLOSE(54.717029611, result1.value().lat().v(), 100 * 1e-3);
  BOOST_CHECK_CLOSE(-14.56385575, result1.value().lon().v(), 100 * 1e-3);

  // Swap geodesics
  const auto result2{
      calculate_intersection_point(g2, g1, units::si::Metres(1e-3))};
  BOOST_CHECK(result2.has_value());

  BOOST_CHECK_CLOSE(54.717029611, result2.value().lat().v(), 100 * 1e-3);
  BOOST_CHECK_CLOSE(-14.56385575, result2.value().lon().v(), 100 * 1e-3);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_closest_intersection_lengths_baselga_3) {
  // Third example from Baselga paper
  // GeodesicArc<double> line_1(35.0, -92.0, 40.0, 52.0);
  // GeodesicArc<double> line_2(-8.0, 20.0, 49.0, -95.0);
  const LatLong a1(Degrees(35.0), Degrees(-92.0));
  const LatLong b1(Degrees(40.0), Degrees(52.0));
  const LatLong a2(Degrees(-8.0), Degrees(20.0));
  const LatLong b2(Degrees(49.0), Degrees(-95.0));

  const Geodesic<double> g1(a1, b1);
  const Geodesic<double> g2(a2, b2);

  const auto result1{
      calculate_intersection_point(g1, g2, units::si::Metres(1e-3))};
  BOOST_CHECK(result1.has_value());

  BOOST_CHECK_CLOSE(50.47909744, result1.value().lat().v(), 100 * 1e-3);
  BOOST_CHECK_CLOSE(-79.2828016944, result1.value().lon().v(), 100 * 1e-3);

  // Swap geodesics
  const auto result2{
      calculate_intersection_point(g2, g1, units::si::Metres(1e-3))};
  BOOST_CHECK(result2.has_value());

  BOOST_CHECK_CLOSE(50.47909744, result2.value().lat().v(), 100 * 1e-3);
  BOOST_CHECK_CLOSE(-79.2828016944, result2.value().lon().v(), 100 * 1e-3);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_closest_intersection_point_non_intersecting) {
  const LatLong istanbul(Degrees(42.0), Degrees(29.0));
  const LatLong washington(Degrees(39.0), Degrees(-77.0));
  const LatLong reyjavik(Degrees(64.0), Degrees(-22.0));
  const LatLong accra(Degrees(6.0), Degrees(0.0));

  const Geodesic<double> g1(istanbul, accra);
  const Geodesic<double> g2(reyjavik, washington);

  const auto result1{
      calculate_intersection_point(g1, g2, units::si::Metres(1e-3))};
  BOOST_CHECK(!result1.has_value());
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
