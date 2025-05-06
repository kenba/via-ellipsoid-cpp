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
/// @file test_Ellipsoid.cpp
/// @brief Contains unit tests for the via::ellipsoid Ellipsoid class.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid/Ellipsoid.hpp"
#include <boost/test/unit_test.hpp>

using namespace via::ellipsoid;
using namespace via;

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_Ellipsoid)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Ellipsoid_wgs84_double) {
  // WGS84 constructor
  const auto ellipsoid{Ellipsoid<double>::wgs84()};
  BOOST_CHECK_EQUAL(wgs84::A<double>.v(), ellipsoid.a().v());
  BOOST_CHECK_EQUAL(wgs84::F<double>, ellipsoid.f());

  BOOST_CHECK_EQUAL(calculate_minor_axis(wgs84::A<double>, wgs84::F<double>),
                    ellipsoid.b());

  BOOST_CHECK_EQUAL(1 - wgs84::F<double>, ellipsoid.one_minus_f());
  BOOST_CHECK_EQUAL(1 / (1 - wgs84::F<double>), ellipsoid.recip_one_minus_f());
  BOOST_CHECK_EQUAL(calculate_sq_eccentricity(wgs84::F<double>),
                    ellipsoid.e_2());
  BOOST_CHECK_EQUAL(calculate_sq_2nd_eccentricity(wgs84::F<double>),
                    ellipsoid.ep_2());
  BOOST_CHECK_EQUAL(calculate_3rd_flattening(wgs84::F<double>), ellipsoid.n());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_Ellipsoid_wgs84_long_double) {
  // WGS84 constructor
  const auto ellipsoid{Ellipsoid<long double>::wgs84()};
  BOOST_CHECK_EQUAL(wgs84::A<long double>.v(), ellipsoid.a().v());
  BOOST_CHECK_EQUAL(wgs84::F<long double>, ellipsoid.f());

  BOOST_CHECK_EQUAL(
      calculate_minor_axis(wgs84::A<long double>, wgs84::F<long double>),
      ellipsoid.b());

  BOOST_CHECK_EQUAL(1 - wgs84::F<long double>, ellipsoid.one_minus_f());
  BOOST_CHECK_EQUAL(1 / (1 - wgs84::F<long double>),
                    ellipsoid.recip_one_minus_f());
  BOOST_CHECK_EQUAL(calculate_sq_eccentricity(wgs84::F<long double>),
                    ellipsoid.e_2());
  BOOST_CHECK_EQUAL(calculate_sq_2nd_eccentricity(wgs84::F<long double>),
                    ellipsoid.ep_2());
  BOOST_CHECK_EQUAL(calculate_3rd_flattening(wgs84::F<long double>),
                    ellipsoid.n());
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
