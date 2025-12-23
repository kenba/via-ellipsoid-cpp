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
                                   units::si::Metres(0.0), g.ellipsoid());
  // a geodesic from the mid point of g to its end
  const Angle<double> half_arc_length_angle(half_arc_length);
  const GeodesicSegment<double> g2(
      g.arc_beta(half_arc_length_angle), g.arc_longitude(half_arc_length),
      g.arc_azimuth(half_arc_length_angle), half_arc_length,
      units::si::Metres(0.0), g.ellipsoid());

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
      g.arc_beta(half_arc_length_angle), g.arc_longitude(half_arc_length),
      g.azi(), half_arc_length, units::si::Metres(0.0), g.ellipsoid());

  // geodesics are NOT coincident
  const auto [distance3, distance4, angle2, iterations2]{
      calculate_sphere_intersection_distances(g1, g3, precision)};
  BOOST_CHECK_CLOSE(g1.arc_length().v(), distance3.v(), precision.v());
  BOOST_CHECK_SMALL(distance4.v(), 16 * CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(
      (g3.azi() - g1.arc_azimuth(half_arc_length_angle)).to_degrees().v(),
      angle2.to_degrees().v(), precision.v());
  BOOST_CHECK_EQUAL(5, iterations2);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_same_geodesic_no_intersection) {
  using namespace via::trig;
  const LatLong a(Degrees(0.0), Degrees(-4.0));
  const LatLong b(Degrees(0.0), Degrees(0.0));
  const GeodesicSegment<double> g1(a, b);
  const GeodesicSegment<double> g3(b, a);

  const LatLong c(Degrees(0.0), Degrees(0.25));
  const LatLong d(Degrees(0.0), Degrees(4.0));
  const GeodesicSegment<double> g2(c, d);
  const GeodesicSegment<double> g4(d, c);

  // 1m precision in Radians on the unit sphere
  const Radians<double> precision{1.0 / g1.ellipsoid().a().v()};

  const auto [distance1, distance2, angle1, iterations]{
      calculate_sphere_intersection_distances(g1, g2, precision)};
  BOOST_CHECK_CLOSE(deg2rad(4.2643), distance1.v(), 400 * precision.v());
  BOOST_CHECK_EQUAL(0.0, distance2.v());

  const auto [distance1_1, distance2_1, angle1_1, iterations_1]{
      calculate_sphere_intersection_distances(g3, g2, precision)};
  BOOST_CHECK_CLOSE(deg2rad(-0.250841), distance1_1.v(), 400 * precision.v());
  BOOST_CHECK_EQUAL(0.0, distance2_1.v());

  const auto [distance1_2, distance2_2, angle1_2, iterations_2]{
      calculate_sphere_intersection_distances(g1, g4, precision)};
  BOOST_CHECK_CLOSE(deg2rad(4.2643), distance1_2.v(), 400 * precision.v());
  BOOST_CHECK_EQUAL(g4.arc_length().v(), distance2_2.v());

  const auto [distance1_3, distance2_3, angle1_3, iterations_3]{
      calculate_sphere_intersection_distances(g3, g4, precision)};
  BOOST_CHECK_EQUAL(0.0, distance1_3.v());
  BOOST_CHECK_CLOSE(deg2rad(4.013456), distance2_3.v(), 400 * precision.v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_intersection_point_non_wgs84) {
  // Example from Charles Karney email on 31/03/2025
  const Ellipsoid ellipsoid(units::si::Metres(6.4e6), 1.0 / 50.0);
  const LatLong p1(Degrees(-30.0), Degrees(0.0));
  const LatLong p2(Degrees(29.5), Degrees(179.5));
  const GeodesicSegment<double> g1(p1, p2, units::si::Metres(0.0),
                                   Radians(great_circle::MIN_VALUE<double>),
                                   ellipsoid);
  BOOST_CHECK_CLOSE(19847901.117944598, g1.length().v(), CALCULATION_TOLERANCE);

  const LatLong p3(Degrees(1.0), Degrees(90.0));
  const LatLong p4(Degrees(-2.0), Degrees(-95.0));
  const GeodesicSegment<double> g2(p3, p4, units::si::Metres(0.0),
                                   Radians(great_circle::MIN_VALUE<double>),
                                   ellipsoid);
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
BOOST_AUTO_TEST_CASE(
    test_calculate_arc_reference_distances_and_angle_coincident_great_circles) {
  const LatLong latlong_w1(Degrees(0.0), Degrees(-1.0));
  const LatLong latlong_e1(Degrees(0.0), Degrees(1.0));
  const GeodesicSegment<double> g_0(latlong_w1, latlong_e1);

  // 1m precision in Radians on the unit sphere
  const Radians<double> precision{1.0 / g_0.ellipsoid().a().v()};

  // same segments
  const auto result_0 =
      calculate_arc_reference_distances_and_angle(g_0, g_0, precision);
  BOOST_CHECK_EQUAL(Radians(0.0), get<0>(result_0));
  BOOST_CHECK_EQUAL(Radians(0.0), get<1>(result_0));
  BOOST_CHECK_EQUAL(Degrees(0.0), get<2>(result_0).to_degrees());

  // opposite segments and same geodesic paths
  const LatLong latlong_w179(Degrees(0.0), Degrees(-179.0));
  const LatLong latlong_e179(Degrees(0.0), Degrees(179.0));
  const GeodesicSegment<double> g_1(latlong_e179, latlong_w179);
  const auto result_1 =
      calculate_arc_reference_distances_and_angle(g_0, g_1, precision);
  BOOST_CHECK_CLOSE(-trig::PI_2<double>, get<0>(result_1).v(), precision.v());
  BOOST_CHECK_CLOSE(trig::PI_2<double>, get<1>(result_1).v(), precision.v());
  BOOST_CHECK_EQUAL(Degrees(0.0), get<2>(result_1).to_degrees());

  // opposite segments and geodesic paths
  const GeodesicSegment<double> g_2(latlong_w179, latlong_e179);
  const auto result_2 =
      calculate_arc_reference_distances_and_angle(g_0, g_2, precision);
  BOOST_CHECK_CLOSE(-trig::PI_2<double>, get<0>(result_2).v(), precision.v());
  BOOST_CHECK_CLOSE(-trig::PI_2<double>, get<1>(result_2).v(), precision.v());
  BOOST_CHECK_EQUAL(Degrees(180.0), get<2>(result_2).to_degrees());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(
    test_calculate_arc_reference_distances_and_angle_intersecting_great_circles) {
  const LatLong latlong_w1(Degrees(0.0), Degrees(-1.0));
  const LatLong latlong_e1(Degrees(0.0), Degrees(1.0));
  const GeodesicSegment<double> g_0(latlong_w1, latlong_e1);

  const LatLong latlong_s1(Degrees(-1.0), Degrees(0.0));
  const LatLong latlong_n1(Degrees(1.0), Degrees(0.0));
  const GeodesicSegment<double> g_1(latlong_s1, latlong_n1);

  // 1m precision in Radians on the unit sphere
  const Radians<double> precision{1.0 / g_0.ellipsoid().a().v()};

  // intersection, same mid points, acute angle
  const auto result_0 =
      calculate_arc_reference_distances_and_angle(g_0, g_1, precision);
  BOOST_CHECK_SMALL(get<0>(result_0).v(), precision.v());
  BOOST_CHECK_SMALL(get<1>(result_0).v(), precision.v());
  BOOST_CHECK_EQUAL(Degrees(90.0), get<2>(result_0).to_degrees());

  const auto angle{44.81195977064123};

  const LatLong latlong_sw1(Degrees(-1.0), Degrees(-1.0));
  const LatLong latlong_ne1(Degrees(1.0), Degrees(1.0));
  const GeodesicSegment<double> g_2(latlong_sw1, latlong_ne1);
  const auto result_1 =
      calculate_arc_reference_distances_and_angle(g_0, g_2, precision);
  BOOST_CHECK_SMALL(get<0>(result_1).v(), precision.v());
  BOOST_CHECK_SMALL(get<1>(result_1).v(), precision.v());
  BOOST_CHECK_CLOSE(angle, get<2>(result_1).to_degrees().v(), precision.v());

  // intersection, same mid points, obtuse angle
  const GeodesicSegment<double> g_3(latlong_ne1, latlong_sw1);
  const auto result_2 =
      calculate_arc_reference_distances_and_angle(g_0, g_3, precision);
  BOOST_CHECK_SMALL(get<0>(result_2).v(), precision.v());
  BOOST_CHECK_SMALL(get<1>(result_2).v(), precision.v());
  BOOST_CHECK_CLOSE(180.0 - angle, get<2>(result_2).to_degrees().v(),
                    precision.v());

  // intersection, different mid points, acute angle
  const GeodesicSegment<double> g_4(
      Angle<double>(), Angle<double>(),
      g_2.arc_azimuth(Angle(g_2.arc_length().half())),
      Radians(trig::PI_2<double>), units::si::Metres(0.0), g_2.ellipsoid());
  const auto result_3 =
      calculate_arc_reference_distances_and_angle(g_0, g_4, precision);
  BOOST_CHECK_SMALL(get<0>(result_3).v(), precision.v());
  BOOST_CHECK_CLOSE(-trig::PI_4<double>, get<1>(result_3).v(), precision.v());
  BOOST_CHECK_CLOSE(angle, get<2>(result_3).to_degrees().v(), precision.v());

  // intersection, different mid points, obtuse angle
  const GeodesicSegment<double> g_5{g_4.reverse()};
  const auto result_4 =
      calculate_arc_reference_distances_and_angle(g_0, g_5, precision);
  BOOST_CHECK_SMALL(get<0>(result_4).v(), precision.v());
  BOOST_CHECK_CLOSE(trig::PI_4<double>, get<1>(result_4).v(), precision.v());
  BOOST_CHECK_CLOSE(180.0 - angle, get<2>(result_4).to_degrees().v(),
                    precision.v());
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
