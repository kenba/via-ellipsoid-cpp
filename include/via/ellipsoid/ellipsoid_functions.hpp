#pragma once

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
/// @file ellipsoid_functions.hpp
/// @brief Contains the via::ellipsoid template functions.
//////////////////////////////////////////////////////////////////////////////
#include "Metres.hpp"
#include <via/angle.hpp>

namespace via {
namespace ellipsoid {
/// Calculate the Semiminor axis of an ellipsoid.
/// @pre a < 0
/// @pre f != 0
/// @param a the Semimajor axis of an ellipsoid.
/// @param f the flattening ratio.
/// @return Semiminor axis of an ellipsoid.
template <typename T>
  requires std::floating_point<T>
constexpr auto calculate_minor_axis(const Metres<T> a, const T f) noexcept
    -> Metres<T> {
  return Metres<T>(a.v() * (T(1) - f));
}

/// Calculate the square of the Eccentricity of an ellipsoid.
/// @pre f != 0
/// @param f the flattening ratio.
/// @return square of the Eccentricity of an ellipsoid.
template <typename T>
  requires std::floating_point<T>
constexpr auto calculate_sq_eccentricity(const T f) noexcept -> T {
  return f * (2 - f);
}

/// Calculate the square of the second Eccentricity of an ellipsoid.
/// @pre f != 0
/// @param f the flattening ratio.
/// @return square of the Eccentricity of an ellipsoid.
template <typename T>
  requires std::floating_point<T>
constexpr auto calculate_sq_2nd_eccentricity(const T f) noexcept -> T {
  const T one_minus_f{1 - f};
  return calculate_sq_eccentricity(f) / (one_minus_f * one_minus_f);
}

/// Calculate the third flattening of an ellipsoid.
/// @pre f != 0
/// @param f the flattening ratio.
/// @return third flattening of an ellipsoid.
template <typename T>
  requires std::floating_point<T>
constexpr auto calculate_3rd_flattening(const T f) noexcept -> T {
  return f / (2 - f);
}

/// Function to calculate epsilon, the variable used in series expansions,
/// derived from Clairaut's constant. Note: epsilon is positive and small.
/// @see CFF Karney, Geodesics on an ellipsoid of revolution: Eqs 22 & 46,
/// https://arxiv.org/pdf/1102.1215.pdf.
/// @pre -1 <= clairaut <= 1
/// @param clairaut Clairaut's constant.
/// @param ep_2 the square of the second Eccentricity of the ellipsoid,
/// default square of WGS 84 second Eccentricity.
/// @return epsilon
template <typename T>
  requires std::floating_point<T>
constexpr auto calculate_epsilon(const trig::UnitNegRange<T> clairaut,
                                 const T ep_2) noexcept -> T {
  // Clairaut's constant is sin alpha0; sq_cos_alpha0 is 1 - clairaut^2
  const T sq_cos_alpha0{(1 - clairaut.v()) * (1 + clairaut.v())};
  const T k2{ep_2 * sq_cos_alpha0}; // square of Karney equation 22
  const T sqrt_k2_1{std::sqrt(1 + k2) + 1};
  return k2 / (sqrt_k2_1 * sqrt_k2_1); // Karney equation 46
}

/// Function to convert an geodetic sine Latitude to a parametric sine
/// Latitude on the auxiliary sphere.
/// @pre -1.0 <= sin_lat <= 1.0
/// @post -1.0 <= result <= 1.0
/// @param sin_lat the sine of the geodetic Latitude.
/// @param e_2 the square of the Eccentricity of the ellipsoid.
/// @return the sine of the parametric Latitude on the auxiliary sphere.
template <typename T>
  requires std::floating_point<T>
constexpr auto
calculate_parametric_sin_latitude(const trig::UnitNegRange<T> sin_lat,
                                  const T e_2) noexcept
    -> trig::UnitNegRange<T> {
  return trig::UnitNegRange<T>::clamp(
      sin_lat.v() *
      std::sqrt((1 - e_2) / (1 - e_2 * sin_lat.v() * sin_lat.v())));
}

/// Function to convert a parametric sine Latitude on the auxiliary sphere
/// to the geodetic sine Latitude.
/// @pre -1.0 <= sin_lat <= 1.0
/// @post -1.0 <= result <= 1.0
/// @param sin_lat the sine of the parametric Latitude.
/// @param e_2 the square of the Eccentricity of the ellipsoid.
/// @return the sine of the geodetic Latitude.
template <typename T>
  requires std::floating_point<T>
constexpr auto
calculate_geodetic_sin_latitude(const trig::UnitNegRange<T> sin_lat,
                                const T e_2) noexcept -> trig::UnitNegRange<T> {
  return trig::UnitNegRange<T>::clamp(
      sin_lat.v() / std::sqrt(1 - e_2 * (1 - sin_lat.v()) * (1 + sin_lat.v())));
}

/// Function to convert a `geodetic` Latitude to a `parametric` Latitude on the
/// auxiliary sphere.
/// @param lat the the geodetic Latitude in an Angle.
/// @param one_minus_f one minus the flattening of the ellipsoid.
/// @return the parametric Latitude on the auxiliary sphere in an Angle.
template <typename T>
  requires std::floating_point<T>
constexpr auto calculate_parametric_latitude(const Angle<T> &lat,
                                             const T one_minus_f) noexcept
    -> Angle<T> {
  return Angle<T>(one_minus_f * lat.sin().v(), lat.cos().v());
}

/// Function to convert a `parametric` Latitude on the auxiliary sphere to a
/// `geodetic` Latitude.
/// @param lat the the geodetic Latitude in an Angle.
/// @param one_minus_f one minus the flattening of the ellipsoid.
/// @return the geodetic Latitude in an Angle.
template <typename T>
  requires std::floating_point<T>
constexpr auto calculate_geodetic_latitude(const Angle<T> &lat,
                                           const T one_minus_f) noexcept
    -> Angle<T> {
  return Angle<T>(lat.sin().v() / one_minus_f, lat.cos().v());
}

} // namespace ellipsoid
} // namespace via
