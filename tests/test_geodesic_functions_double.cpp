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
/// @file test_geodesic_functions_double.cpp
/// @brief Contains unit tests for the via::ellipsoid geodesic_functions.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid/GeodesicSegment.hpp"
#include "via/ellipsoid/wgs84.hpp"
#include <boost/test/unit_test.hpp>
#include <via/ellipsoid/Ellipsoid.hpp>

using namespace via::ellipsoid;
using namespace via;

namespace {
const auto EPSILON(std::numeric_limits<double>::epsilon());
const auto CALCULATION_TOLERANCE(100 * EPSILON);
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_geodesic_functions)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_astroid_double) {
  // test zero cases
  BOOST_CHECK_EQUAL(0.0, calculate_astroid(0.0, 0.0));
  BOOST_CHECK_EQUAL(0.0, calculate_astroid(1.0, 0.0));

  // 0.0, 0.0 to 0.5, 179.5
  BOOST_CHECK_CLOSE(
      0.91583665308532092,
      calculate_astroid(-0.82852367684428574, -0.82576675584253256),
      CALCULATION_TOLERANCE);

  // 0.0, 0.0 to 1.0, 179.0
  BOOST_CHECK_CLOSE(1.9858096632693705,
                    calculate_astroid(-1.6572357126833825, -1.6518470456464789),
                    CALCULATION_TOLERANCE);

  // -30.0, 0.0 to 30.0, 179.0
  BOOST_CHECK_CLOSE(0.91211900939748036,
                    calculate_astroid(-1.9121190093974805, 0.0),
                    CALCULATION_TOLERANCE);

  // -30.0, 0.0 to 30.5, 179.5
  BOOST_CHECK_CLOSE(
      1.2324261949931818,
      calculate_astroid(-0.96091919533424308, -1.1124132048023443),
      CALCULATION_TOLERANCE);

  const auto X_THRESHOLD{1000 * std::sqrt(EPSILON)};
  const auto Y_TOLERANCE{200 * EPSILON};

  BOOST_CHECK_CLOSE(1.4901161193852098e-05,
                    calculate_astroid(-1.0 - X_THRESHOLD, -Y_TOLERANCE),
                    CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_end_azimuth_double) {
  const auto angle_90{Angle<double>::from_y_x(1.0, 0.0)};
  const Angle<double> angle_50(Degrees<double>(50.0));
  const auto angle_45{Angle<double>::from_y_x(1.0, 1.0)};
  const Angle<double> angle_20(Degrees<double>(20.0));

  auto result = calculate_end_azimuth(angle_20, angle_50, angle_20);
  BOOST_CHECK_CLOSE(30.0, result.to_degrees().v(), CALCULATION_TOLERANCE);

  result = calculate_end_azimuth(angle_50, angle_20, angle_20);
  BOOST_CHECK_CLOSE(13.530064432438888, result.to_degrees().v(),
                    CALCULATION_TOLERANCE);

  result = calculate_end_azimuth(-angle_50, angle_50, angle_20);
  BOOST_CHECK_EQUAL(20.0, result.to_degrees().v());

  result = calculate_end_azimuth(angle_45, angle_45, angle_90);
  BOOST_CHECK_EQUAL(90.0, result.to_degrees().v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_delta_omega12_double) {
  constexpr double WGS84_EP_2{calculate_sq_2nd_eccentricity(wgs84::F<double>)};

  // 0.0, 0.0 to 30.0, 90.0
  const auto clairaut_30_90{Angle(Degrees(60.0)).sin()};
  const auto lam12_30_90{delta_omega12(
      clairaut_30_90, calculate_epsilon(clairaut_30_90, WGS84_EP_2),
      Radians(trig::PI_2<double>), Angle<double>(), Angle(Degrees(90.0)),
      Ellipsoid<double>::wgs84())};
  BOOST_CHECK_CLOSE(0.0045600360192803542, lam12_30_90.v(),
                    CALCULATION_TOLERANCE);

  // 0.0, 0.0 to 45.0, 90.0
  const auto clairaut_45_90{Angle(Degrees(45.0)).sin()};
  const auto lam12_45_90{delta_omega12(
      clairaut_45_90, calculate_epsilon(clairaut_45_90, WGS84_EP_2),
      Radians(trig::PI_2<double>), Angle<double>(), Angle(Degrees(90.0)),
      Ellipsoid<double>::wgs84())};
  BOOST_CHECK_CLOSE(0.0037224722989948442, lam12_45_90.v(),
                    2 * CALCULATION_TOLERANCE);

  // 0.0, 0.0 to 60.0, 90.0
  const auto clairaut_60_90{Angle(Degrees(30.0)).sin()};
  const auto lam12_60_90{delta_omega12(
      clairaut_60_90, calculate_epsilon(clairaut_60_90, WGS84_EP_2),
      Radians(trig::PI_2<double>), Angle<double>(), Angle(Degrees(90.0)),
      Ellipsoid<double>::wgs84())};
  BOOST_CHECK_CLOSE(0.0026316334829412581, lam12_60_90.v(),
                    CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_azimuths_arc_length_meridian_double) {
  const LatLong latlon1(Degrees(-70.0), Degrees(40.0));
  const LatLong latlon2(Degrees(80.0), Degrees(40.0));

  // Northbound geodesic along a meridian
  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
  BOOST_CHECK_EQUAL(0.0, azimuth1.to_degrees().v());
  BOOST_CHECK_CLOSE(2.6163378712682306, arc_length1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(0.0, end_azimuth1.to_degrees().v());
  BOOST_CHECK_EQUAL(0, iter1);

  // Southbound geodesic along a meridian
  const auto [azimuth2, arc_length2, end_azimuth2,
              iter2]{calculate_azimuths_arc_length(latlon2, latlon1)};
  BOOST_CHECK_EQUAL(180.0, azimuth2.to_degrees().v());
  BOOST_CHECK_CLOSE(2.6163378712682306, arc_length2.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(180.0, end_azimuth2.to_degrees().v());
  BOOST_CHECK_EQUAL(0, iter2);

  // Northbound geodesic past the North pole
  const LatLong latlon3(Degrees(80.0), Degrees(-140.0));
  const auto [azimuth3, arc_length3, end_azimuth3,
              iter3]{calculate_azimuths_arc_length(latlon2, latlon3)};
  BOOST_CHECK_EQUAL(0.0, azimuth3.to_degrees().v());
  BOOST_CHECK_CLOSE(0.3502163200513691, arc_length3.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(180.0, end_azimuth3.to_degrees().v());
  BOOST_CHECK_EQUAL(0, iter3);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_azimuths_arc_length_equator_double) {
  const LatLong latlon1(Degrees(0.0), Degrees(-40.0));
  const LatLong latlon2(Degrees(0.0), Degrees(50.0));

  // Eastbound geodesic along the equator
  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
  BOOST_CHECK_EQUAL(90.0, azimuth1.to_degrees().v());
  BOOST_CHECK_CLOSE(1.5760806267286946, arc_length1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(90.0, end_azimuth1.to_degrees().v());
  BOOST_CHECK_EQUAL(0, iter1);

  // Westbound geodesic along the equator
  const auto [azimuth2, arc_length2, end_azimuth2,
              iter2]{calculate_azimuths_arc_length(latlon2, latlon1)};
  BOOST_CHECK_EQUAL(-90.0, azimuth2.to_degrees().v());
  BOOST_CHECK_CLOSE(1.5760806267286946, arc_length2.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(-90.0, end_azimuth2.to_degrees().v());
  BOOST_CHECK_EQUAL(0, iter2);

  // Long Eastbound geodesic along the equator
  const LatLong latlon3(Degrees(0.0), Degrees(135.0));
  const auto [azimuth3, arc_length3, end_azimuth3,
              iter3]{calculate_azimuths_arc_length(latlon1, latlon3)};
  BOOST_CHECK_EQUAL(90.0, azimuth3.to_degrees().v());
  BOOST_CHECK_CLOSE(3.0646012186391296, arc_length3.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(90.0, end_azimuth3.to_degrees().v());
  BOOST_CHECK_EQUAL(0, iter3);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(
    test_calculate_azimuths_arc_length_equator_nearly_antipodal_double) {
  const LatLong latlon1(Degrees(0.0), Degrees(0.0));
  const LatLong latlon2(Degrees(0.0), Degrees(179.5));

  // Eastbound geodesic along the equator
  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
#ifdef _MSC_VER
  BOOST_CHECK_CLOSE(55.966495140158635, azimuth1.to_degrees().v(),
                    2 * CALCULATION_TOLERANCE);
#else
  BOOST_CHECK_EQUAL(55.966495140158635, azimuth1.to_degrees().v());
#endif
  BOOST_CHECK_EQUAL(trig::PI<double>, arc_length1.v());
#ifdef _MSC_VER
  BOOST_CHECK_CLOSE(180.0 - 55.966495140158635, end_azimuth1.to_degrees().v(),
                    2 * CALCULATION_TOLERANCE);
#else
  BOOST_CHECK_EQUAL(180.0 - 55.966495140158635, end_azimuth1.to_degrees().v());
#endif
  BOOST_CHECK_EQUAL(3, iter1);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_arc_length_aziumth_normal_01) {
  const LatLong latlon1(Degrees(-40.0), Degrees(70.0));
  const LatLong latlon2(Degrees(30.0), Degrees(0.0));

  // North West bound, straddle Equator
  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
  BOOST_CHECK_CLOSE(-55.00473169905792, azimuth1.to_degrees().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(1.6656790467428877, arc_length1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(-46.47061016713593, end_azimuth1.to_degrees().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(3, iter1);

  const auto beta1{Ellipsoid<double>::wgs84().calculate_parametric_latitude(
      Angle(Degrees(-40.0)))};
  const GeodesicSegment<double> g(beta1, Angle(Degrees(70.0)), azimuth1,
                                  arc_length1);
  const auto latlon3{g.arc_lat_long(arc_length1, Angle(arc_length1))};
  BOOST_CHECK_CLOSE(latlon2.lat().v(), latlon3.lat().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_SMALL(latlon3.lon().v(), CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_arc_length_aziumth_normal_02) {
  const LatLong latlon1(Degrees(30.0), Degrees(70.0));
  const LatLong latlon2(Degrees(-40.0), Degrees(0.0));

  // South West bound, straddle Equator
  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
  BOOST_CHECK_EQUAL(-133.52938983286407, azimuth1.to_degrees().v());
  BOOST_CHECK_CLOSE(1.6656790467428877, arc_length1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(-124.99526830094207, end_azimuth1.to_degrees().v());
  BOOST_CHECK_EQUAL(3, iter1);

  const auto beta1{Ellipsoid<double>::wgs84().calculate_parametric_latitude(
      Angle(Degrees(30.0)))};
  const GeodesicSegment<double> g(beta1, Angle(Degrees(70.0)), azimuth1,
                                  arc_length1);
  const auto latlon3{g.arc_lat_long(arc_length1, Angle(arc_length1))};
  BOOST_CHECK_CLOSE(latlon2.lat().v(), latlon3.lat().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_SMALL(latlon3.lon().v(), CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_arc_length_aziumth_normal_03) {
  const LatLong latlon1(Degrees(30.0), Degrees(0.0));
  const LatLong latlon2(Degrees(-40.0), Degrees(70.0));

  // South East bound, straddle Equator
  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
  BOOST_CHECK_EQUAL(133.52938983286407, azimuth1.to_degrees().v());
  BOOST_CHECK_CLOSE(1.6656790467428877, arc_length1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(124.99526830094207, end_azimuth1.to_degrees().v());
  BOOST_CHECK_EQUAL(3, iter1);

  const auto beta1{Ellipsoid<double>::wgs84().calculate_parametric_latitude(
      Angle(Degrees(30.0)))};
  const GeodesicSegment<double> g(beta1, Angle(Degrees(0.0)), azimuth1,
                                  arc_length1);
  const auto latlon3{g.arc_lat_long(arc_length1, Angle(arc_length1))};
  BOOST_CHECK_CLOSE(latlon2.lat().v(), latlon3.lat().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(latlon2.lon().v(), latlon3.lon().v(),
                    CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_arc_length_aziumth_normal_04) {
  const LatLong latlon1(Degrees(-40.0), Degrees(0.0));
  const LatLong latlon2(Degrees(30.0), Degrees(70.0));

  // North East bound, straddle Equator
  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
  BOOST_CHECK_CLOSE(55.00473169905792, azimuth1.to_degrees().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(1.6656790467428874, arc_length1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(46.470610167135938, end_azimuth1.to_degrees().v());
  BOOST_CHECK_EQUAL(3, iter1);

  const auto beta1{Ellipsoid<double>::wgs84().calculate_parametric_latitude(
      Angle(Degrees(-40.0)))};
  const GeodesicSegment<double> g(beta1, Angle(Degrees(0.0)), azimuth1,
                                  arc_length1);
  const auto latlon3{g.arc_lat_long(arc_length1, Angle(arc_length1))};
  BOOST_CHECK_CLOSE(latlon2.lat().v(), latlon3.lat().v(),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(latlon2.lon().v(), latlon3.lon().v(),
                    CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_arc_length_aziumth_normal_05) {
  const LatLong latlon1(Degrees(0.0), Degrees(0.0));
  const LatLong latlon2(Degrees(0.5), Degrees(179.98));

  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
  BOOST_CHECK_CLOSE(1.0420381519981656,
                    azimuth1.to_degrees().v(), // 1.0420381519981552
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(3.132893826005981, arc_length1.v(), CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(178.9579224301469, end_azimuth1.to_degrees().v());
  BOOST_CHECK_EQUAL(3, iter1);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_geodesic_arc_length_aziumth_normal_06) {
  // GeodTest.dat line 460107
  // 89.985810803742 0 90.033923043742 -89.985810803761488692
  // 179.999716989078075251 89.966210133068275597 20003931.4528694
  // 179.999999966908046132 .0036837809003 -47969483155.576793
  const LatLong latlon1(Degrees(89.985810803742), Degrees(0.0));
  const LatLong latlon2(Degrees(-89.985810803761488692),
                        Degrees(179.999716989078075251));

  const auto [azimuth1, arc_length1, end_azimuth1,
              iter1]{calculate_azimuths_arc_length(latlon1, latlon2)};
  BOOST_CHECK_CLOSE(90.033923043742, azimuth1.to_degrees().v(), 1.7e-5);
  BOOST_CHECK_CLOSE(89.966210133068275597, end_azimuth1.to_degrees().v(),
                    1.7e-5);

  const auto beta1{Ellipsoid<double>::wgs84().calculate_parametric_latitude(
      Angle(Degrees(89.985810803742)))};
  const auto distance{convert_radians_to_metres(beta1, azimuth1, arc_length1)};
  BOOST_CHECK_CLOSE(20003931.4528694, distance.v(), CALCULATION_TOLERANCE);
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
