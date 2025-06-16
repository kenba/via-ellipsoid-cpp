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
/// @file test_GeodesicSegment_double.cpp
/// @brief Contains unit tests for the via::ellipsoid GeodesicSegment class.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid/GeodesicSegment.hpp"
#include <boost/test/unit_test.hpp>
#include <via/angle/trig.hpp>

using namespace via::ellipsoid;
using namespace via;

namespace {
const auto EPSILON(std::numeric_limits<double>::epsilon());
const auto CALCULATION_TOLERANCE(100 * EPSILON);
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_Geodesic_double)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Geodesic_direct_constructors) {
  const units::si::Metres<double> length{9'000'000.0};
  const Radians arc_length{trig::PI_2<double>};

  // Ensure that two Geodesics can fit on a cache line.
  BOOST_CHECK_EQUAL(128u, sizeof(GeodesicSegment<double>));

  // Latitude, longitude, azimuth, "direct" constructor.
  const GeodesicSegment<double> geodesic0(Angle(Degrees(45.0)),
                                          Angle(Degrees(0.0)),
                                          Angle(Degrees(90.0)), Radians(0.0));
  BOOST_CHECK_EQUAL(
      trig::PI_2<double>,
      geodesic0.arc_azimuth(Angle(Radians(0.0))).to_radians().v());

  // LatLong, azimuth, arc_length constructor.
  const LatLong a(Degrees(45.0), Degrees(45.0));

  // Increase azimuth around compass from due South to due North
  for (auto i{-179}; i <= 180; ++i) {
    const Angle<double> azimuth{Degrees<double>(i)};
    const GeodesicSegment<double> geodesic1(a, azimuth, length);
    BOOST_CHECK(geodesic1.is_valid());

    const Angle<double> azi0{geodesic1.azimuth(units::si::Metres(0.0))};
    BOOST_CHECK_CLOSE(trig::deg2rad<double>(i), azi0.to_radians().v(),
                      2 * CALCULATION_TOLERANCE);

    const auto len0{geodesic1.length()};
    BOOST_CHECK_CLOSE(length.v(), len0.v(), CALCULATION_TOLERANCE);

    const GeodesicSegment<double> geodesic2(a, azimuth, arc_length);
    BOOST_CHECK(geodesic2.is_valid());

    const Angle<double> azi2{geodesic2.arc_azimuth(Angle(Radians(0.0)))};
    BOOST_CHECK_CLOSE(trig::deg2rad<double>(i), azi2.to_radians().v(),
                      2 * CALCULATION_TOLERANCE);

    const auto len2{geodesic2.arc_length()};
    BOOST_CHECK_CLOSE(arc_length.v(), len2.v(), CALCULATION_TOLERANCE);
  }
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Geodesic_between_positions) {
  const LatLong istanbul(Degrees(42.0), Degrees(29.0));
  const LatLong washington(Degrees(39.0), Degrees(-77.0));

  const GeodesicSegment<double> g1(istanbul, washington);
  BOOST_CHECK(g1.is_valid());

  const auto end_azimuth{g1.azimuth(g1.length())};
  BOOST_CHECK_CLOSE(-132.2646607116376, end_azimuth.to_degrees().v(),
                    CALCULATION_TOLERANCE);

  // test start position
  BOOST_CHECK_CLOSE(
      istanbul.lat().v(),
      g1.ellipsoid().calculate_geodetic_latitude(g1.beta()).to_degrees().v(),
      CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(istanbul.lon().v(), g1.lon().to_degrees().v(),
                    CALCULATION_TOLERANCE);

  const auto lat_long{g1.lat_long(units::si::Metres(0.0))};
  BOOST_CHECK_CLOSE(istanbul.lat().v(), lat_long.lat().v(),
                    2 * CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(istanbul.lon().v(), lat_long.lon().v(),
                    CALCULATION_TOLERANCE);

  // test start point
  const auto start_point{g1.arc_point(Radians(0.0))};
  BOOST_CHECK_CLOSE(
      istanbul.lat().v(),
      g1.ellipsoid()
          .calculate_geodetic_latitude(vector::latitude(start_point))
          .to_degrees()
          .v(),
      CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(istanbul.lon().v(),
                    vector::longitude(start_point).to_degrees().v(),
                    CALCULATION_TOLERANCE);

  const auto [atd, xtd, iterations]{
      g1.calculate_atd_and_xtd(istanbul, units::si::Metres(1e-3))};
  BOOST_CHECK_EQUAL(0.0, atd.v());
  BOOST_CHECK_EQUAL(0.0, xtd.v());

  auto distamce{g1.shortest_distance(istanbul, units::si::Metres(1e-3))};
  BOOST_CHECK_EQUAL(0.0, distamce.v());

  // test end position
  const auto arc_length{g1.arc_length()};
  BOOST_CHECK_CLOSE(1.309412846249522, arc_length.v(), CALCULATION_TOLERANCE);

  const auto length{g1.length()};
  BOOST_CHECK_CLOSE(8339863.136005359, length.v(), CALCULATION_TOLERANCE);

  BOOST_CHECK_CLOSE(washington.lat().v(), g1.latitude(length).to_degrees().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(washington.lon().v(), g1.longitude(length).to_degrees().v(),
                    CALCULATION_TOLERANCE);

  // test mid position
  const units::si::Metres<double> half_length{0.5 * g1.length().v()};
  const auto mid_position{g1.lat_long(half_length)};

  BOOST_CHECK_CLOSE(54.86379153725445, mid_position.lat().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(-25.694568908316413, mid_position.lon().v(),
                    4 * CALCULATION_TOLERANCE);

  const auto mid_length{g1.metres_to_radians(half_length)};
  BOOST_CHECK_CLOSE(0.654673165141749, mid_length.v(), CALCULATION_TOLERANCE);
  const auto mid_point{g1.mid_point()};
  const auto mid_beta{vector::latitude(mid_point)};
  const auto mid_lat{g1.ellipsoid().calculate_geodetic_latitude(mid_beta)};
  BOOST_CHECK_CLOSE(mid_position.lat().v(), mid_lat.to_degrees().v(),
                    CALCULATION_TOLERANCE);
  const auto mid_lon{vector::longitude(mid_point)};
  BOOST_CHECK_CLOSE(mid_position.lon().v(), mid_lon.to_degrees().v(),
                    CALCULATION_TOLERANCE);

  const Radians<double> precision(1e-3 / wgs84::A<double>.v());
  const auto [atd2, xtd2, iterations2]{
      g1.calculate_arc_atd_and_xtd(mid_beta, mid_lon, precision)};
  BOOST_CHECK_CLOSE(mid_length.v(), atd2.v(), 100 * precision.v());
  BOOST_CHECK_SMALL(xtd2.v(), 100 * precision.v());

  const auto [atd3, xtd3, iterations3]{
      g1.calculate_atd_and_xtd(mid_position, units::si::Metres(1e-3))};
  BOOST_CHECK_CLOSE(half_length.v(), atd3.v(), 100 * precision.v());
  BOOST_CHECK_SMALL(xtd3.v(), 100 * precision.v());

  distamce = g1.shortest_distance(mid_position, units::si::Metres(1e-3));
  BOOST_CHECK_EQUAL(0.0, distamce.v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_meridonal_Geodesic) {
  const LatLong a(Degrees(45.0), Degrees(0.0));
  const LatLong b(Degrees(45.0), Degrees(180.0));
  const GeodesicSegment<double> g1(a, b);
  BOOST_CHECK(g1.is_valid());
  const auto [_point, pole0]{g1.arc_point_and_pole(Radians(0.0))};

  // Calculate the azimuth at the North pole
  const Radians mid_length(0.5 * g1.arc_length().v());
  const auto azimuth{g1.arc_azimuth(Angle(mid_length))};
  BOOST_CHECK_EQUAL(180.0, azimuth.to_degrees().v());

  // Calculate the point and great circle pole at the North pole
  const auto [point1, pole1]{g1.arc_point_and_pole(mid_length)};
  BOOST_CHECK_EQUAL(pole0, pole1);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Geodesic_90n_0n_0e) {
  const LatLong a(Degrees(90.0), Degrees(0.0));
  const LatLong b(Degrees(0.0), Degrees(0.0));
  const GeodesicSegment<double> g1(a, b);
  BOOST_CHECK(g1.is_valid());

  BOOST_CHECK_CLOSE(trig::PI_2<double>, g1.arc_length().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(180.0, g1.azi().to_degrees().v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Geodesic_90s_0n_50e) {
  const LatLong a(Degrees(-90.0), Degrees(0.0));
  const LatLong b(Degrees(0.0), Degrees(50.0));
  const GeodesicSegment<double> g1(a, b);
  BOOST_CHECK(g1.is_valid());

  BOOST_CHECK_CLOSE(trig::PI_2<double>, g1.arc_length().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(0.0, g1.azi().to_degrees().v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_atd_and_xtd) {
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

  // Test delta_azimuth at interception, should be PI/2
  const auto azimuth_1{g1.azimuth(atd)};
  const GeodesicSegment<double> g2(position, reyjavik);
  const auto azimuth_2{g2.azi()};
  const auto delta_azimuth{azimuth_2 - azimuth_1};
  BOOST_CHECK_CLOSE(trig::PI_2<double>, delta_azimuth.to_radians().v(),
                    200 * CALCULATION_TOLERANCE);

  // opposite geodesic
  const GeodesicSegment<double> g3(washington, istanbul);
  const auto [atd2, xtd2, iterations2]{
      g3.calculate_atd_and_xtd(reyjavik, units::si::Metres(1e-3))};
  BOOST_CHECK_CLOSE(g1.length().v() - 3928788.572, atd2.v(), 100 * 1e-3);
  BOOST_CHECK_CLOSE(1010585.9988368, xtd2.v(), 100 * 1e-3);

  // const LatLong reyjavik(Degrees(64.0), Degrees(-22.0));
  auto distamce{g1.shortest_distance(reyjavik, units::si::Metres(1e-3))};
  BOOST_CHECK_CLOSE(1010585.998836817, distamce.v(), 1e-3);

  const LatLong accra(Degrees(6.0), Degrees(0.0));
  distamce = g1.shortest_distance(accra, units::si::Metres(1e-3));
  BOOST_CHECK_CLOSE(4891211.398445355, distamce.v(), 1e-3);

  const LatLong chicago(Degrees(42.0), Degrees(-88.0));
  distamce = g1.shortest_distance(chicago, units::si::Metres(1e-3));
  BOOST_CHECK_CLOSE(989277.1859906457, distamce.v(), 1e-3);

  const LatLong singapore(Degrees(1.0), Degrees(104.0));
  distamce = g1.shortest_distance(singapore, units::si::Metres(1e-3));
  BOOST_CHECK_CLOSE(8699538.22763653, distamce.v(), 1e-3);
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
