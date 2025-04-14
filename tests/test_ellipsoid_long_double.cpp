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
/// @file test_ellipsoid_long_double.cpp
/// @brief Contains unit tests for the via::ellipsoid Geodesic class.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace via::ellipsoid;
using namespace via;

namespace {
const auto EPSILON(std::numeric_limits<long double>::epsilon());
const auto CALCULATION_TOLERANCE(100 * EPSILON);
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_ellipsoid_long_double)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_intersection_point_non_wgs84) {
  // Example from Charles Karney email on 01/04/2025
  const LatLong p1(Degrees(4.0L), Degrees(20.0L));
  const Angle a1(Degrees(-56.0L));
  const Geodesic<long double> g1(p1, a1, Radians(trig::PI<long double> - 0.1L));

  const LatLong p2(Degrees(-30.0L), Degrees(-40.0L));
  const Angle a2(Degrees(80.0L));
  const Geodesic<long double> g2(p2, a2, Radians(trig::PI<long double> - 0.1L));

  const auto result{calculate_aux_intersection_distances(
      g1, g2, Radians(great_circle::MIN_VALUE<long double>))};
  const auto d1{std::get<0>(result)};
  const auto d2{std::get<1>(result)};
#ifdef OUTPUT_ITERATIONS
  const auto iterations{std::get<2>(result)};
  std::cout << "Long double precision (Radians): "
            << great_circle::MIN_VALUE<long double> << std::endl;
  std::cout << "Long double iterations: " << iterations << std::endl;
#endif

  const auto pd1{g1.aux_lat_long(d1)};
  const auto pd2{g2.aux_lat_long(d2)};
// Disable tests since Visual Studio cannot calculate long doubles
#ifndef _MSC_VER
  BOOST_CHECK_CLOSE(-1.44147956008236583L, pd1.lat().v(),
                    16 * CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(27.97257917717199337L, pd1.lon().v(),
                    2 * CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(pd1.lat().v(), pd2.lat().v(), 18 * CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(pd1.lon().v(), pd2.lon().v(), CALCULATION_TOLERANCE);
#endif
#ifdef OUTPUT_ITERATIONS
  // calculate number of iterations for 1mm precision
  const auto result_2{calculate_aux_intersection_distances(
      g1, g2, Radians<long double>(1e-3 / g1.ellipsoid().a().v()))};
  const auto iterations_2{std::get<2>(result_2)};
  std::cout << "1mm precision (Radians): " << 1e-3 / g1.ellipsoid().a().v()
            << std::endl;
  std::cout << "1mm iterations: " << iterations_2 << std::endl;
#endif
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
