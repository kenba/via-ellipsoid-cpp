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
/// @file test_ellipsoid_functions_double.cpp
/// @brief Contains unit tests for the via::ellipsoid template functions.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid/ellipsoid_functions.hpp"
#include "via/ellipsoid/wgs84.hpp"
#include <boost/test/unit_test.hpp>

using namespace via::ellipsoid;
using namespace via;

namespace {
const auto EPSILON(std::numeric_limits<double>::epsilon());
const auto CALCULATION_TOLERANCE(100 * EPSILON);
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_ellipsoid_functions)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_epsilon) {
  const auto wgs84_ep2{calculate_sq_2nd_eccentricity(wgs84::F<double>)};
  BOOST_CHECK_CLOSE(0.0016792203863837047,
                    calculate_epsilon(trig::UnitNegRange(0.0), wgs84_ep2),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(0.0015745990877544997,
                    calculate_epsilon(trig::UnitNegRange(0.25), wgs84_ep2),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(0.0012604720416530619,
                    calculate_epsilon(trig::UnitNegRange(0.5), wgs84_ep2),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_CLOSE(0.0007360477262034019,
                    calculate_epsilon(trig::UnitNegRange(0.75), wgs84_ep2),
                    CALCULATION_TOLERANCE);
  BOOST_CHECK_EQUAL(0.0, calculate_epsilon(trig::UnitNegRange(1.0), wgs84_ep2));
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_calculate_parametric_and_geodetic_latitude) {
  const auto one_minus_f{1.0 - wgs84::F<double>};

  for (auto latitude(-90); latitude <= 90; ++latitude) {
    const auto lat{Angle(Degrees<double>(latitude))};
    const auto parametric_lat{calculate_parametric_latitude(lat, one_minus_f)};
    const auto result(calculate_geodetic_latitude(parametric_lat, one_minus_f));

    BOOST_CHECK_CLOSE(lat.to_radians().v(), result.to_radians().v(),
                      2 * CALCULATION_TOLERANCE);
    // The largest differences occur at 3 and 50 degrees.
  }
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
