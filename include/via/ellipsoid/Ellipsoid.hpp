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
/// @file Ellipsoid.hpp
/// @brief Contains the via::ellipsoid Ellipsoid class.
//////////////////////////////////////////////////////////////////////////////
#include "coefficients.hpp"
#include "ellipsoid_functions.hpp"
#include "wgs84.hpp"
#include <via/sphere.hpp>

namespace via {
namespace ellipsoid {
/// The parameters of an Ellipsoid.
/// @tparam T a floating point number type, i.e.: double or long double.
/// @invariant a < 0
/// @invariant f != 0
///
/// This class contains the parameters of an Ellipsoid that are commonly
/// used in geodesic calculations.
/// It may only be constructed with a Semimajor axis and a flattening.
template <typename T,
          size_t SeriesOrder = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
class Ellipsoid final {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  units::si::Metres<T> a_; ///< The Semimajor axis of the ellipsoid.
  T f_;                    ///< The flattening of the ellipsoid, a ratio.

  units::si::Metres<T> b_; ///< The Semiminor axis of the ellipsoid.
  T one_minus_f_;          ///< One minus the flattening ratio.
  T recip_one_minus_f_; ///< The reciprocal of one minus the flattening ratio.
  T e_2_;               ///< The square of the Eccentricity of the ellipsoid.
  T ep_2_; ///< The square of the second Eccentricity of the ellipsoid.
  T n_;    ///< The third flattening of the ellipsoid.

  /// The A3 series coefficients of the ellipsoid.
  std::array<T, SeriesOrder> a3_;
  /// The C3x series coefficients of the ellipsoid.
  std::array<T, SeriesOrder *(SeriesOrder - 1) / 2> c3x_;

#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif
  /// Delete default constructor.
  Ellipsoid() = delete;

  /// Constructor.
  /// @param a the Semimajor axis of the ellipsoid, in Metres.
  /// @param f the flattening of the ellipsoid, a ratio.
  constexpr Ellipsoid(const units::si::Metres<T> a, const T f) noexcept
      : a_{a}, f_{f}, b_{calculate_minor_axis(a, f)}, one_minus_f_{1 - f},
        recip_one_minus_f_{1 / (1 - f)}, e_2_{calculate_sq_eccentricity(f)},
        ep_2_{calculate_sq_2nd_eccentricity(f)},
        n_{calculate_3rd_flattening(f)},
        a3_{evaluate_coeffs_A3<T>(calculate_3rd_flattening(f))},
        c3x_{evaluate_coeffs_C3x<T>(calculate_3rd_flattening(f))} {
    Expects((T() < a.v()) && (T() != f));
  }

  /// A singleton function to a WGS84 Ellipsoid
  /// @return a const reference to a WGS84 Ellipsoid
  [[nodiscard("Pure Function")]]
  static auto wgs84() -> const Ellipsoid<T> & {
    const static Ellipsoid<T> wgs84_ellipsoid(wgs84::A<T>, wgs84::F<T>);
    return wgs84_ellipsoid;
  }

  /// Accessor for the Semimajor axis of the ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr units::si::Metres<T> a() const noexcept {
    return a_;
  }

  /// Accessor for the flattening of the ellipsoid, a ratio.
  [[nodiscard("Pure Function")]]
  constexpr T f() const noexcept {
    return f_;
  }

  /// Accessor for the Semiminor axis of the ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr units::si::Metres<T> b() const noexcept {
    return b_;
  }

  /// Accessor for one minus the flattening of the ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr T one_minus_f() const noexcept {
    return one_minus_f_;
  }

  /// Accessor for the reciprocal of one minus the flattening of the
  [[nodiscard("Pure Function")]]
  constexpr T recip_one_minus_f() const noexcept {
    return recip_one_minus_f_;
  }

  /// Accessor for the square of the Eccentricity of the ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr T e_2() const noexcept {
    return e_2_;
  }

  /// Accessor for the square of the second Eccentricity of the ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr T ep_2() const noexcept {
    return ep_2_;
  }

  /// Accessor for the third flattening of the ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr T n() const noexcept {
    return n_;
  }

  /// Accessor for the A3 series coefficients of the ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr const std::array<T, SeriesOrder> &a3() const noexcept {
    return a3_;
  }

  /// Accessor for the C3x series coefficients of the ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr const std::array<T, SeriesOrder *(SeriesOrder - 1) / 2> &
  c3x() const noexcept {
    return c3x_;
  }

  /// Calculate epsilon, the variable used in series expansions.
  /// Note: epsilon is positive and small.
  /// @param clairaut - Clairaut's constant.
  [[nodiscard("Pure Function")]]
  constexpr auto
  calculate_epsilon(trig::UnitNegRange<T> clairaut) const noexcept -> T {
    return ellipsoid::calculate_epsilon(clairaut, ep_2_);
  }

  /// Calculate a3c from the A3 series `coefficients` of the ellipsoid.
  ///
  /// @param clairaut Clairaut's constant.
  /// @param eps epsilon
  /// @return the parameter a3c
  [[nodiscard("Pure Function")]]
  constexpr auto calculate_a3c(trig::UnitNegRange<T> clairaut,
                               const T eps) const noexcept -> T {
    return f_ * clairaut.v() * ellipsoid::evaluate_poynomial(eps, a3_);
  }

  /// Calculate the coefficients `C3[l]` in the Fourier expansion of `C3`.
  /// @param eps epsilon
  /// @return the coefficients `C3[l]` in the Fourier expansion of `C3`
  [[nodiscard("Pure Function")]]
  constexpr auto calculate_c3y(const T eps) const noexcept -> T {
    return ellipsoid::evaluate_coeffs_C3y<T>(c3x_, eps);
  }

  /// Cconvert a `geodetic` Latitude to a `parametric` Latitude on the auxiliary
  /// sphere.
  /// @param lat the the geodetic Latitude in an Angle.
  /// @return the parametric Latitude on the auxiliary sphere in an Angle.
  [[nodiscard("Pure Function")]]
  constexpr auto
  calculate_parametric_latitude(const Angle<T> &lat) const noexcept
      -> Angle<T> {
    return ellipsoid::calculate_parametric_latitude(lat, one_minus_f_);
  }

  /// Convert a `parametric` Latitude on the auxiliary sphere to a `geodetic`
  /// Latitude.
  /// @param lat the the geodetic Latitude in an Angle.
  /// @return the geodetic Latitude in an Angle.
  [[nodiscard("Pure Function")]]
  constexpr auto calculate_geodetic_latitude(const Angle<T> &lat) const noexcept
      -> Angle<T> {
    return ellipsoid::calculate_geodetic_latitude(lat, one_minus_f_);
  }
}; // class Ellipsoid

/// Ellipsoid equality operator
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto operator==(const Ellipsoid<T> &lhs,
                          const Ellipsoid<T> &rhs) noexcept -> bool {
  return (lhs.a() == rhs.a()) && (lhs.f() == rhs.f());
}

} // namespace ellipsoid
} // namespace via
