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
/// @file test_ellipsoid_double.cpp
/// @brief Contains unit tests for the via::ellipsoid GeodesicSegment class.
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
  const GeodesicSegment<double> g1(istanbul, washington);
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

  const GeodesicSegment<double> g1(a1, b1);
  const GeodesicSegment<double> g2(a2, b2);

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

  const GeodesicSegment<double> g1(istanbul, washington);
  const GeodesicSegment<double> g2(reyjavik, accra);

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

  const GeodesicSegment<double> g1(a1, b1);
  const GeodesicSegment<double> g2(a2, b2);

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

  const GeodesicSegment<double> g1(istanbul, accra);
  const GeodesicSegment<double> g2(reyjavik, washington);

  const auto result1{
      calculate_intersection_point(g1, g2, units::si::Metres(1e-3))};
  BOOST_CHECK(!result1.has_value());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_intersection_same_geodesic_split) {
  const LatLong a(Degrees(1.0), Degrees(0.0));
  const LatLong b(Degrees(-0.998286322222), Degrees(179.296674991667));
  const GeodesicSegment<double> g(a, b);

  // Split g into two geodesics
  const units::si::Metres<double> half_length{g.length().v() / 2.0};
  const auto half_arc_length{g.metres_to_radians(half_length)};

  // a geodesic from the start of g to its mid point
  const GeodesicSegment<double> g1(g.beta(), g.lon(), g.azi(), half_arc_length,
                                   g.ellipsoid());
  // a geodesic from the mid point of g to its end
  const Angle<double> half_arc_length_angle(half_arc_length);
  const GeodesicSegment<double> g2(
      g.arc_beta(half_arc_length_angle),
      g.arc_longitude(half_arc_length, half_arc_length_angle),
      g.arc_azimuth(half_arc_length_angle), half_arc_length, g.ellipsoid());

  // 1mm precision in Radians on the auxiliary sphere
  const Radians<double> precision{1e-3 / g.ellipsoid().a().v()};

  // geodesics are coincident
  const auto [distance1, distance2, angle1, iterations]{
      calculate_sphere_intersection_distances(g1, g2, precision)};
  BOOST_CHECK_CLOSE(g1.arc_length().v(), distance1.v(),
                    2 * CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(0.0, distance2.v());
  BOOST_CHECK_EQUAL(0.0, angle1.to_degrees().v());
  BOOST_CHECK_EQUAL(0, iterations);

  // a geodesic from the mid point of g to another point
  const GeodesicSegment<double> g3(
      g.arc_beta(half_arc_length_angle),
      g.arc_longitude(half_arc_length, half_arc_length_angle), g.azi(),
      half_arc_length, g.ellipsoid());

  // geodesics are NOT coincident
  const auto [distance3, distance4, angle2, iterations2]{
      calculate_sphere_intersection_distances(g1, g3, precision)};
  BOOST_CHECK_CLOSE(g1.arc_length().v(), distance3.v(),
                    2 * CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(0.0, distance4.v());
  BOOST_CHECK_EQUAL(0.0, angle2.to_degrees().v());
  BOOST_CHECK_EQUAL(0, iterations2);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_intersection_point_non_wgs84) {
  // Example from Charles Karney email on 31/03/2025
  const Ellipsoid ellipsoid(units::si::Metres(6.4e6), 1.0 / 50.0);
  const LatLong p1(Degrees(-30.0), Degrees(0.0));
  const LatLong p2(Degrees(29.5), Degrees(179.5));
  const GeodesicSegment<double> g1(
      p1, p2, Radians(great_circle::MIN_VALUE<double>), ellipsoid);
  BOOST_CHECK_CLOSE(19847901.117944598, g1.length().v(), CALCULATION_TOLERANCE);

  const LatLong p3(Degrees(1.0), Degrees(90.0));
  const LatLong p4(Degrees(-2.0), Degrees(-95.0));
  const GeodesicSegment<double> g2(
      p3, p4, Radians(great_circle::MIN_VALUE<double>), ellipsoid);
  BOOST_CHECK_CLOSE(19518466.043349568, g2.length().v(), CALCULATION_TOLERANCE);

  const auto result{
      calculate_intersection_point(g1, g2, units::si::Metres(1e-6))};
  BOOST_CHECK_CLOSE(-28.099944988083493, result->lat().v(),
                    24 * CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(172.27633238700983, result->lon().v(),
                    CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
