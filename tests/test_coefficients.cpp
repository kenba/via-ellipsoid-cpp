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
/// @file test_coefficients.hpp
/// @brief Contains unit tests for the via::ellipsoid coefficients and
/// series expansion functions.
//////////////////////////////////////////////////////////////////////////////
#include "via/ellipsoid/coefficients.hpp"
#include "via/ellipsoid/ellipsoid_functions.hpp"
#include <boost/geometry/util/series_expansion.hpp>
#include <boost/test/unit_test.hpp>

using namespace via;
using namespace via::ellipsoid;

namespace {
constexpr double WGS84_F{1.0 / 298.257223563};
constexpr long double WGS84_F_L{1.0L / 298.257223563L};

// The third flattening of the WGS84 ellipsoid.
constexpr double N{calculate_3rd_flattening(WGS84_F)};
constexpr long double N_L{calculate_3rd_flattening(WGS84_F_L)};

// Corresponds to the Clairaut's constant for 45 degrees azimuth
constexpr double EPS{calculate_sq_2nd_eccentricity(WGS84_F) / 2};
constexpr long double EPS_L{calculate_sq_2nd_eccentricity(WGS84_F_L) / 2};
} // namespace

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(Test_coefficients)

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_A1) {
  auto boost_A1_6(boost::geometry::series_expansion::evaluate_A1<6>(EPS));
  auto via_A1_6(evaluate_A1(EPS));
  BOOST_CHECK_EQUAL(boost_A1_6, via_A1_6);

  auto boost_A1_8(boost::geometry::series_expansion::evaluate_A1<8>(EPS_L));
  auto via_A1_8(evaluate_A1<long double, 8>(EPS_L));
  BOOST_CHECK_EQUAL(boost_A1_8, via_A1_8);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_A2) {
  auto boost_A2_6(boost::geometry::series_expansion::evaluate_A2<6>(EPS));
  auto via_A2_6(evaluate_A2(EPS));
  BOOST_CHECK_EQUAL(boost_A2_6, via_A2_6);

  auto boost_A2_8(boost::geometry::series_expansion::evaluate_A2<8>(EPS_L));
  auto via_A2_8(evaluate_A2<long double, 8>(EPS_L));
  BOOST_CHECK_EQUAL(boost_A2_8, via_A2_8);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_coeffs_A3) {
  boost::array<double, 6> boost_A3_6;
  boost::geometry::series_expansion::evaluate_coeffs_A3(boost_A3_6, N);
  auto via_A3_6(evaluate_coeffs_A3<double>(N));
  for (auto i(0u); i < 6u; ++i)
    BOOST_CHECK_EQUAL(boost_A3_6[i], via_A3_6[i]);

  boost::array<long double, 8> boost_A3_8;
  boost::geometry::series_expansion::evaluate_coeffs_A3(boost_A3_8, N_L);
  auto via_A3_8(evaluate_coeffs_A3<long double, 8>(N_L));
  for (auto i(0u); i < 8u; ++i)
    BOOST_CHECK_EQUAL(boost_A3_8[i], via_A3_8[i]);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_coeffs_C1) {
  boost::array<double, 7> boost_C1_6;
  boost::geometry::series_expansion::evaluate_coeffs_C1(boost_C1_6, EPS);
  auto via_C1_6(evaluate_coeffs_C1<double>(EPS));
  for (auto i(1u); i < 7u; ++i)
    BOOST_CHECK_EQUAL(boost_C1_6[i], via_C1_6[i]);

  boost::array<long double, 9> boost_C1_8;
  boost::geometry::series_expansion::evaluate_coeffs_C1(boost_C1_8, EPS_L);
  auto via_C1_8(evaluate_coeffs_C1<long double, 8>(EPS_L));
  for (auto i(1u); i < 9u; ++i)
    BOOST_CHECK_EQUAL(boost_C1_8[i], via_C1_8[i]);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_coeffs_C1p) {
  boost::array<double, 6> boost_C1p_6;
  boost::geometry::series_expansion::evaluate_coeffs_C1p(boost_C1p_6, EPS);
  auto via_C1p_6(evaluate_coeffs_C1p<double>(EPS));
  for (auto i(1u); i < 6u; ++i)
    BOOST_CHECK_EQUAL(boost_C1p_6[i], via_C1p_6[i]);

  boost::array<long double, 8> boost_C1p_8;
  boost::geometry::series_expansion::evaluate_coeffs_C1p(boost_C1p_8, EPS_L);
  auto via_C1p_8(evaluate_coeffs_C1p<long double, 8>(EPS_L));
  for (auto i(1u); i < 8u; ++i)
    BOOST_CHECK_EQUAL(boost_C1p_8[i], via_C1p_8[i]);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_coeffs_C2) {
  boost::array<double, 7> boost_C2_6;
  boost::geometry::series_expansion::evaluate_coeffs_C2(boost_C2_6, EPS);
  auto via_C2_6(evaluate_coeffs_C2<double>(EPS));
  for (auto i(1u); i < 7u; ++i)
    BOOST_CHECK_EQUAL(boost_C2_6[i], via_C2_6[i]);

  boost::array<long double, 9> boost_C2_8;
  boost::geometry::series_expansion::evaluate_coeffs_C2(boost_C2_8, EPS_L);
  auto via_C2_8(evaluate_coeffs_C2<long double, 8>(EPS_L));
  for (auto i(1u); i < 9u; ++i)
    BOOST_CHECK_EQUAL(boost_C2_8[i], via_C2_8[i]);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_coeffs_C3x) {
  boost::array<double, 15> boost_C3x_6;
  boost::geometry::series_expansion::evaluate_coeffs_C3x<6>(boost_C3x_6, N);
  auto via_C3x_6(evaluate_coeffs_C3x<double>(N));
  for (auto i(0u); i < 15u; ++i)
    BOOST_CHECK_EQUAL(boost_C3x_6[i], via_C3x_6[i]);

  boost::array<long double, 28> boost_C3x_8;
  boost::geometry::series_expansion::evaluate_coeffs_C3x<8>(boost_C3x_8, N_L);
  auto via_C3x_8(evaluate_coeffs_C3x<long double, 8>(N_L));
  for (auto i(0u); i < 28u; ++i)
    BOOST_CHECK_EQUAL(boost_C3x_8[i], via_C3x_8[i]);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_coeffs_C3y) {
  boost::array<double, 15> boost_C3x_6;
  boost::geometry::series_expansion::evaluate_coeffs_C3x<6>(boost_C3x_6, N);
  boost::array<double, 6> boost_C3_6;
  boost::geometry::series_expansion::evaluate_coeffs_C3(boost_C3_6, boost_C3x_6,
                                                        EPS);

  auto via_C3x_6(evaluate_coeffs_C3x<double>(N));
  auto via_C3_6(evaluate_coeffs_C3y<double>(via_C3x_6, EPS));
  for (auto i(1u); i < 6u; ++i)
    BOOST_CHECK_EQUAL(boost_C3_6[i], via_C3_6[i]);

  boost::array<long double, 28> boost_C3x_8;
  boost::geometry::series_expansion::evaluate_coeffs_C3x<8>(boost_C3x_8, N_L);
  boost::array<long double, 8> boost_C3_8;
  boost::geometry::series_expansion::evaluate_coeffs_C3(boost_C3_8, boost_C3x_8,
                                                        EPS_L);

  auto via_C3x_8(evaluate_coeffs_C3x<long double, 8>(N_L));
  auto via_C3_8(evaluate_coeffs_C3y<long double, 8>(via_C3x_8, EPS_L));
  for (auto i(1u); i < 8u; ++i)
    BOOST_CHECK_EQUAL(boost_C3_8[i], via_C3_8[i]);
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_sin_cos_series_C3) {
  const Angle<double> sigma(Radians(0.1 * M_PI));

  boost::array<double, 15> boost_C3x_6;
  boost::geometry::series_expansion::evaluate_coeffs_C3x<6>(boost_C3x_6, N);
  boost::array<double, 6> boost_C3_6;
  boost::geometry::series_expansion::evaluate_coeffs_C3(boost_C3_6, boost_C3x_6,
                                                        EPS);
  auto boost_C3_6_B31(boost::geometry::series_expansion::sin_cos_series(
      sigma.sin().v(), sigma.cos().v(), boost_C3_6));

  auto via_C3s_6(evaluate_coeffs_C3x(N));
  auto via_C3_6(evaluate_coeffs_C3y(via_C3s_6, EPS));
  auto via_C3_6_B31(sin_cos_series(sigma, via_C3_6));
  BOOST_CHECK_EQUAL(boost_C3_6_B31, via_C3_6_B31.v());

  boost::array<long double, 28> boost_C3x_8;
  boost::geometry::series_expansion::evaluate_coeffs_C3x<8>(boost_C3x_8, N_L);
  boost::array<long double, 8> boost_C3_8;
  boost::geometry::series_expansion::evaluate_coeffs_C3(boost_C3_8, boost_C3x_8,
                                                        EPS_L);
  auto boost_C3_8_B31(boost::geometry::series_expansion::sin_cos_series(
      sigma.sin().v(), sigma.cos().v(), boost_C3_8));

  auto via_C3x_8(evaluate_coeffs_C3x<long double, 8>(N_L));
  auto via_C3_8(evaluate_coeffs_C3y<long double, 8>(via_C3x_8, EPS_L));
  auto via_C3_8_B31(sin_cos_series(sigma, via_C3_8));
  BOOST_CHECK_EQUAL(boost_C3_8_B31, via_C3_8_B31.v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_sin_cos_series_C3_2) {
  const Angle<double> sigma(Radians(0.1 * M_PI));

  boost::array<double, 15> boost_C3x_6;
  boost::geometry::series_expansion::evaluate_coeffs_C3x<6>(boost_C3x_6, N);
  boost::array<double, 6> boost_C3_6;
  boost::geometry::series_expansion::evaluate_coeffs_C3(boost_C3_6, boost_C3x_6,
                                                        EPS);
  auto boost_C3_6_B31(boost::geometry::series_expansion::sin_cos_series(
      sigma.sin().v(), sigma.cos().v(), boost_C3_6));

  auto via_C3s_6(evaluate_coeffs_C3x(N));
  auto via_C3_6(evaluate_coeffs_C3y(via_C3s_6, EPS));
  auto via_C3_6_B31(sin_cos_series(sigma, via_C3_6));
  BOOST_CHECK_EQUAL(boost_C3_6_B31, via_C3_6_B31.v());

  boost::array<long double, 28> boost_C3x_8;
  boost::geometry::series_expansion::evaluate_coeffs_C3x<8>(boost_C3x_8, N_L);
  boost::array<long double, 8> boost_C3_8;
  boost::geometry::series_expansion::evaluate_coeffs_C3(boost_C3_8, boost_C3x_8,
                                                        EPS_L);
  auto boost_C3_8_B31(boost::geometry::series_expansion::sin_cos_series(
      sigma.sin().v(), sigma.cos().v(), boost_C3_8));

  auto via_C3x_8(evaluate_coeffs_C3x<long double, 8>(N_L));
  auto via_C3_8(evaluate_coeffs_C3y<long double, 8>(via_C3x_8, EPS_L));
  auto via_C3_8_B31(sin_cos_series(sigma, via_C3_8));
  BOOST_CHECK_EQUAL(boost_C3_8_B31, via_C3_8_B31.v());
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(test_evaluate_poynomial_A3) {
  boost::array<double, 6> boost_A3_6;
  boost::geometry::series_expansion::evaluate_coeffs_A3(boost_A3_6, N);

  auto boost_A3_6_EPS(boost::geometry::math::horner_evaluate(
      EPS, boost_A3_6.cbegin(), boost_A3_6.cend()));
  auto via_A3_6(evaluate_coeffs_A3<double>(N));
  auto via_A3_6_EPS(evaluate_poynomial(EPS, via_A3_6));
  BOOST_CHECK_EQUAL(boost_A3_6_EPS, via_A3_6_EPS);

  boost::array<long double, 8> boost_A3_8;
  boost::geometry::series_expansion::evaluate_coeffs_A3(boost_A3_8, N);
  auto boost_A3_8_EPS(boost::geometry::math::horner_evaluate(
      EPS, boost_A3_8.cbegin(), boost_A3_8.cend()));

  auto via_A3_8(evaluate_coeffs_A3<long double, 8>(N));
  auto via_A3_8_EPS(evaluate_poynomial(EPS, via_A3_8));
  BOOST_CHECK_EQUAL(boost_A3_8_EPS, via_A3_8_EPS);
}
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////////////////////////////////////////
