#pragma once

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
/// @file geodesic_functions.hpp
/// @brief Contains the via::ellipsoid geodesic template functions.
//////////////////////////////////////////////////////////////////////////////
#include "Ellipsoid.hpp"
#include <cmath>
#include <limits>
#include <via/angle/trig.hpp>
#ifdef OUTPUT_GEOD_ITERATOR_STEPS
#include <iomanip>
#include <iostream>
#endif

namespace via {
namespace ellipsoid {

/// Estimate omega12 by solving the astroid problem.
/// Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
/// @param x, y astroid parameters, see Karney section 7.
/// @return k the solution of the astroid.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_astroid(const T x, const T y) -> T {
  const T p{x * x};
  const T q{y * y};
  const T r{(p + q - 1) / T(6)};

  // y = 0 with |x| <= 1
  // for y small, positive root is k = abs(y)/sqrt(1-x^2)
  if ((q == T()) && (r <= T())) {
    return T();
  }

  const T S{p * q / T(4)};
  const T r2{r * r};
  const T r3{r * r2};
  T u{r};

  // The discriminant of the quadratic equation for T3.
  // This is zero on the evolute curve p^(1/3)+q^(1/3) = 1
  const T discriminant{S * (S + 2 * r3)};
  if (discriminant >= T()) {
    T t3{S + r3};
    // Pick the sign on the sqrt to maximize abs(T3), to minimise loss
    // of precision due to cancellation.
    t3 += std::copysign(std::sqrt(discriminant), t3);
    const T t{std::cbrt(t3)};
    u += t + ((t == T()) ? T() : r2 / t);
  } else {
    // discriminant < T()
    // T is complex, but the way u is defined the result is real.
    const T angle{std::atan2(std::sqrt(-discriminant), -(S + r3))};
    // There are three possible cube roots.  We choose the root which
    // avoids cancellation.  Note: discriminant < 0 implies that r < 0.
    u += 2 * r * std::cos(angle / T(3));
  }

  const T v{std::sqrt(u * u + q)};         // guaranteed positive
  const T uv{u < 0 ? q / (v - u) : u + v}; // u+v, guaranteed positive
  const T w{(uv - q) / (2 * v)};           // positive?

  // Rearrange expression for k to avoid loss of accuracy due to subtraction.
  // Division by 0 not possible because uv > 0, w >= 0.
  return uv / (std::sqrt(uv + w * w) + w); // guaranteed positive
}

/// Estimate the initial azimuth for a nearly antipodal arc.
/// It calculates and solves the astroid problem.
///
/// @pre abs_lambda12 >= 0
///
/// @param beta1, beta2 parametric latitudes on the auxiliary sphere.
/// @param abs_lambda12 longitude difference between start and finish points.
/// @param ellipsoid the `Ellipsoid`.
/// @return an estimate of the initial azimuth on the auxiliary sphere.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto estimate_antipodal_initial_azimuth(const Angle<T> beta1,
                                                  const Angle<T> beta2,
                                                  const Angle<T> abs_lambda12,
                                                  const Ellipsoid<T> &ellipsoid)
    -> Angle<T> {
#if defined(_MSC_VER) || defined(__clang__)
  const T X_THRESHOLD{1000 * std::sqrt(std::numeric_limits<T>::epsilon())};
#else
  constexpr T X_THRESHOLD{1000 * std::sqrt(std::numeric_limits<T>::epsilon())};
#endif
  constexpr T Y_TOLERANCE{200 * std::numeric_limits<T>::epsilon()};

  Expects(abs_lambda12.sin().v() >= T());

  // Calculate the integration parameter for geodesic
  const auto clairaut{beta1.cos()}; // Note: assumes sin_alpha_1 = 1
  const T eps{ellipsoid.calculate_epsilon(clairaut)};
  const T a3f{ellipsoid.calculate_a3f(eps)};

  const T lamscale{T(ellipsoid.f()) * beta1.cos().v() * a3f * trig::PI<T>};
  const T betscale{lamscale * beta1.cos().v()};

  // Solve astroid problem
  const T x{abs_lambda12.opposite().to_radians().v() / lamscale};
  const T y{(beta1 + beta2).sin().v() / betscale};

  // Test x and y params
  if ((y > -Y_TOLERANCE) && (x > T(-1) - X_THRESHOLD)) {
    const auto sin_alpha{trig::UnitNegRange<T>::clamp(-x)};
    return Angle<T>(sin_alpha, trig::swap_sin_cos(sin_alpha)).negate_cos();
  } else {
    const T k{calculate_astroid(x, y)};
    const Radians<T> omg12a{lamscale * (-x * k / (1 + k))};

    Angle<T> omega12(omg12a);
    return great_circle::calculate_gc_azimuth(beta1, beta2,
                                              omega12.negate_cos());
  }
}

/// Calculate the azimuth at parametric latitude beta2 from the azimuth
/// at parametric latitude beta1.
///
/// @pre alpha1.sin() > 0
///
/// @param beta1, beta2 parametric latitudes on the auxiliary sphere.
/// @param alpha1 initial azimuth.
/// @return azimuth at beta2.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_end_azimuth(const Angle<T> beta1, const Angle<T> beta2,
                                     const Angle<T> alpha1) -> Angle<T> {
  Expects(alpha1.sin().v() > T());

  const trig::UnitNegRange<T> clairaut{alpha1.sin().v() * beta1.cos().v()};

  const auto sin_alpha2{
      beta2.cos() == beta1.cos()
          ? alpha1.sin()
          : trig::UnitNegRange<T>::clamp(clairaut.v() / beta2.cos().v())};

  // Karney's method to calculate the cosine of the end azimuth
  trig::UnitNegRange<T> cos_alpha2{alpha1.cos().abs()};
  if ((beta2.cos() != beta1.cos()) || (beta2.sin().abs() != -beta1.sin())) {
    const T temp1{alpha1.cos().v() * beta1.cos().v()};
    const T temp2{(beta1.cos() < beta1.sin().abs())
                      ? (beta2.cos().v() - beta1.cos().v()) *
                            (beta1.cos().v() + beta2.cos().v())
                      : (beta1.sin().v() - beta2.sin().v()) *
                            (beta1.sin().v() + beta2.sin().v())};
    const T temp3{temp1 * temp1 + temp2};
    const T temp4{(temp3 > T()) ? std::sqrt(temp3) / beta2.cos().v() : T()};
    cos_alpha2 = trig::UnitNegRange<T>::clamp(temp4);
  }

  return Angle<T>{sin_alpha2, cos_alpha2};
}

/// Calculate the longitude difference between the auxiliary sphere and
/// ellipsoid.
///
/// @pre 0 < sigma12
///
/// @param clairaut Clairaut's constant.
/// @param eps the integration parameter.
/// @param sigma12 great circle distance on the auxiliary sphere.
/// @param sigma1, sigma2 great circle distances on the auxiliary sphere
/// to the points from the Northward Equator crossing.
/// @param ellipsoid the `Ellipsoid`.
/// @return the longitude difference between the auxiliary sphere and
/// ellipsoid in radians.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto delta_omega12(const trig::UnitNegRange<T> clairaut, const T eps,
                             const Radians<T> sigma12, const Angle<T> sigma1,
                             const Angle<T> sigma2,
                             const Ellipsoid<T> &ellipsoid) -> Radians<T> {
  const auto c3{evaluate_coeffs_C3y<T>(ellipsoid.c3x(), eps)};
  const Radians<T> b31{sin_cos_series(sigma1, c3)};
  const Radians<T> b32{sin_cos_series(sigma2, c3)};

  const T a3c{ellipsoid.calculate_a3c(clairaut, eps)};
  return Radians<T>(a3c * (sigma12 + (b32 - b31)).v());
}

/// Estimate the initial azimuth on the auxiliary sphere for normal arcs,
/// i.e. NOT nearly antipodal points.
///
/// @pre abs_lambda12 >= 0
///
/// @param beta1, beta2 parametric latitudes on the auxiliary sphere.
/// @param abs_lambda12 longitude difference between the points.
/// @param ellipsoid the `Ellipsoid`.
/// @return an estimate of the initial azimuth.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto estimate_initial_azimuth(const Angle<T> beta1,
                                        const Angle<T> beta2,
                                        const Angle<T> abs_lambda12,
                                        const Ellipsoid<T> &ellipsoid)
    -> Angle<T> {
  Expects(abs_lambda12.sin().v() >= T());

  // Calculate azimuths at the arc ends
  Angle<T> alpha1{
      great_circle::calculate_gc_azimuth(beta1, beta2, abs_lambda12)};
  Angle<T> alpha2{calculate_end_azimuth(beta1, beta2, alpha1)};

  // Calculate arc lengths from the nearest great circle equator crossing
  const auto sigma1{
      Angle<T>::from_y_x(beta1.sin().v(), beta1.cos().v() * alpha1.cos().v())};
  const auto sigma2{
      Angle<T>::from_y_x(beta2.sin().v(), beta2.cos().v() * alpha2.cos().v())};
  const Radians<T> sigma12{sigma2.to_radians() - sigma1.to_radians()};

  // Calculate Clairaut's constant
  const trig::UnitNegRange<T> clairaut(alpha1.sin().v() * beta1.cos().v());
  const T eps{ellipsoid.calculate_epsilon(clairaut)};

  // Estimate the difference between sphere and geodesic longitude
  const Angle<T> domg12{
      delta_omega12(clairaut, eps, sigma12, sigma1, sigma2, ellipsoid)};
  const auto omega12{abs_lambda12 + domg12};

  // Recalculate the azimuth using the estimated geodesic longitude
  return great_circle::calculate_gc_azimuth(beta1, beta2, omega12);
}

/// Calculate: m12b = (reduced length)/_b
/// @pre sigma12 > 0
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_reduced_length(const T eps, const Radians<T> sigma12,
                                        const Angle<T> sigma1, const T dn1,
                                        const Angle<T> sigma2, const T dn2)
    -> T {
  Expects(sigma12.v() > T());

  const T A1{evaluate_A1<T>(eps)};
  const T A2{evaluate_A2<T>(eps)};
  const T m0x{A1 - A2};

  const T a1p1 = T(1) + A1;
  const T a2p1 = T(1) + A2;

  const auto Ca{evaluate_coeffs_C1<T>(eps)};
  auto Cb{evaluate_coeffs_C2<T>(eps)};

  // Assume here that Ca.size() >= Cb.size()
  // int size(Cb.size() - 1);
  for (auto i(1u); i < Cb.size(); ++i) {
    Cb[i] = a1p1 * Ca[i] - a2p1 * Cb[i];
  }

  const T J12{m0x * sigma12.v() + (sin_cos_series(sigma2, Cb).v() -
                                   sin_cos_series(sigma1, Cb).v())};
  return dn2 * (sigma1.cos().v() * sigma2.sin().v()) -
         dn1 * (sigma1.sin().v() * sigma2.cos().v()) -
         sigma1.cos().v() * sigma2.cos().v() * J12;
}

/// Find the initial aziumth and arc length on the auxiliary sphere given
/// estimates of the initial aziumth: alpha1, and arc length: sigma12.
/// It uses Newton's method to solve:
///   f(alp1) = lambda12(alp1) - lam12 = 0
///
/// @pre beta2 >= beta1
/// @pre abs_lambda12.sin() >= 0
/// @pre alpha1.sin() > 0
/// @pre 0 < sigma12 <= π
/// @post alpha1.sin() > 0
/// @post 0 < sigma12 <= π
///
/// @tparam MAX_ITER1 the first iterations maximum threshold.
/// @param beta1, beta2 parametric latitudes on the auxiliary sphere.
/// @param abs_lambda12 the geodesic longitude difference.
/// @param alpha1 the initial azimuth.
/// @param sigma12 the initial arc length.
/// @param tolerance the tolerance to perform the calculation to in radians.
/// @param ellipsoid the `Ellipsoid`.
/// @return the azimuths at the start of the arc on the auxiliary sphere
///  and arc length on the auxiliary sphere
/// at the start of the geodesic segment and the number of iterations required
/// to calculate them.
template <typename T, unsigned MAX_ITER1 = 20>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto find_azimuth_arc_length_newtons_method(const Angle<T> beta1,
                                            const Angle<T> beta2,
                                            const Angle<T> abs_lambda12,
                                            Angle<T> alpha1, Radians<T> sigma12,
                                            const Radians<T> tolerance,
                                            const Ellipsoid<T> &ellipsoid)
    -> std::tuple<Angle<T>, Radians<T>, unsigned> {
  // maximum iteration threshold from GeographicLib
  constexpr unsigned MAX_ITERS{MAX_ITER1 + std::numeric_limits<T>::digits +
                               10u};
  Expects(beta2 >= beta1);
  Expects(abs_lambda12.sin().v() >= T());
  Expects(alpha1.sin().v() > T());
  Expects((T() < sigma12.v()) && (sigma12.v() <= trig::PI<T>));

  const T dn1{std::sqrt(T(1) + T(ellipsoid.ep_2()) *
                                   (beta1.sin().v() * beta1.sin().v()))};
  const T dn2{std::sqrt(T(1) + T(ellipsoid.ep_2()) *
                                   (beta2.sin().v() * beta2.sin().v()))};

#ifdef OUTPUT_GEOD_ITERATOR_STEPS
  auto prev_v{sigma12.v()};
#endif
  auto iterations{1U};
  for (; iterations <= MAX_ITERS; ++iterations) {
    // Calculate Clairaut's constant
    const trig::UnitNegRange<T> clairaut(alpha1.sin().v() * beta1.cos().v());
    const T eps{ellipsoid.calculate_epsilon(clairaut)};

    // Calculate first longitude (omega1) and distance (sigma1) from the
    // Northbound equator crossing
    const trig::UnitNegRange<T> sin_omega1{clairaut.v() * beta1.sin().v()};
    const auto cos_omega1{beta1.cos().v() * alpha1.cos().v()};
    const auto omega1{Angle<T>::from_y_x(sin_omega1.v(), cos_omega1)};
    const auto sigma1{Angle<T>::from_y_x(beta1.sin().v(), cos_omega1)};

    // Calculate azimuth at the end point
    const Angle<T> alpha2{calculate_end_azimuth(beta1, beta2, alpha1)};

    // Calculate second longitude (omega2) and distance (sigma2) from the
    // Northbound equator crossing
    const trig::UnitNegRange<T> sin_omega2{clairaut.v() * beta2.sin().v()};
    const auto cos_omega2{beta2.cos().v() * alpha2.cos().v()};
    const auto omega2{Angle<T>::from_y_x(sin_omega2.v(), cos_omega2)};
    const auto sigma2{Angle<T>::from_y_x(beta2.sin().v(), cos_omega2)};

    // Calculate great circle length on the auxiliary sphere
    const Angle<T> sc_sigma12{sigma2 - sigma1};
    // clamp to range 0 to Pi
    sigma12 =
        std::signbit(sc_sigma12.sin().v())
            ? (std::signbit(sc_sigma12.cos().v()) ? Radians<T>(trig::PI<T>)
                                                  : Radians<T>())
            : sc_sigma12.to_radians();
    const auto domg12{
        delta_omega12(clairaut, eps, sigma12, sigma1, sigma2, ellipsoid)};

    // Calculate Longitude difference on the auxiliary sphere
    auto omega12 = omega2 - omega1;
    // clamp to range 0 to Pi
    if (std::signbit(omega12.sin().v())) {
      omega12 =
          (std::signbit(omega12.cos().v()))
              ? Angle<T>(trig::UnitNegRange(T()), trig::UnitNegRange(T(-1)))
              : Angle<T>();
    }
    const auto eta{(omega12 - abs_lambda12).to_radians()};

    // Difference between differences
    const T v{eta.v() - domg12.v()};

    // Test within tolerance
    if (std::abs(v) <= tolerance.v()) {
      break;
    }

    // Calculate the denominator for Newton's method
    const T dv{(alpha2.cos().abs().v() == 0.0)
                   ? -2 * ellipsoid.one_minus_f() * dn1 / beta1.sin().v()
                   : ellipsoid.one_minus_f() *
                         calculate_reduced_length(eps, sigma12, sigma1, dn1,
                                                  sigma2, dn2) /
                         (alpha2.cos().v() * beta2.cos().v())};

    // Calculate the change in initial azimuth and test within MIN_VALUE<T>
    const T dalpha1{std::clamp<T>(-v / dv, -1, 1)};

#ifdef OUTPUT_GEOD_ITERATOR_STEPS
    const bool converging{std::abs(v) < prev_v};
    prev_v = std::abs(v);

    std::cout << std::boolalpha;
    std::cout << iterations << ',' << converging << ',' << dalpha1 << ',' << v
              << ',' << dv << std::endl;
#endif

    // Adjust the azimuth by dalpha1 for the next iteration
    alpha1 += Angle<T>(Radians<T>{dalpha1});
  }

  Ensures(alpha1.sin().v() > T());
  Ensures((T() < sigma12.v()) && (sigma12.v() <= trig::PI<T>));

  return {alpha1, sigma12, iterations};
}

/// Find the aziumth and arc length on the auxiliary sphere.
/// It adjusts the latitudes and longitude difference so that the aziumth of
/// the geodesic lies between 0° and 90°.
/// It calls find_azimuth_arc_length_newtons_method and then changes the
/// resulting azimuth to match the orienation of the geodesic.
///
/// @pre beta_a.cos() >= 0
/// @pre beta_b.cos() >= 0
/// @post 0 < arc_length <= π
/// @pre tolerance >= epsilon
/// @post 0 < sigma12 <= π
///
/// @param beta_a, beta_b the geodetic latitudes of the start and finish points.
/// @param lambda12 the geodesic longitude difference in radians.
/// @param arc_length the auxiliary sphere great circle length in radians
/// @param tolerance the tolerance to perform the calculation to in Radians.
/// @param estimate_azimuth a flag to call `estimate_initial_azimuth` before
/// performing Newton's method.
/// @param ellipsoid the `Ellipsoid`.
/// @return the azimuth at the start of the geodesic segment, the great circle
/// arc length on the auxiliary sphere, the azimuth at the end of
/// geodesic segment and the number of iterations required to calculate them.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto find_azimuths_and_arc_length(const Angle<T> beta_a, const Angle<T> beta_b,
                                  const Angle<T> lambda12,
                                  const Radians<T> arc_length,
                                  const Radians<T> tolerance,
                                  const bool estimate_azimuth,
                                  const Ellipsoid<T> &ellipsoid)
    -> std::tuple<Angle<T>, Radians<T>, Angle<T>, unsigned> {
  const T ANTIPODAL_ARC_THRESHOLD{trig::PI<T> * ellipsoid.one_minus_f()};

  Expects(beta_a.cos().v() >= T());
  Expects(beta_b.cos().v() >= T());
  Expects((T() < arc_length.v()) && (arc_length.v() <= trig::PI<T>));
  Expects(tolerance.v() >= std::numeric_limits<T>::epsilon());

  // Start at the latitude furthest from the Equator
  const Angle<T> abs_beta_a{beta_a.abs()};
  const Angle<T> abs_beta_b{beta_b.abs()};
  // Note: the algorithm is very sensitive to starting with the largest latitude
  // so ensure furthest point is used by comparing sines AND cosines
  const bool swap_latitudes{(abs_beta_a.sin() < abs_beta_b.sin()) ||
                            (abs_beta_a.cos() > abs_beta_b.cos())};
  Angle<T> beta1{swap_latitudes ? beta_b : beta_a};
  Angle<T> beta2{swap_latitudes ? beta_a : beta_b};

  // Start South of the Equator
  // Note: sets negate_latitude on the Equator to favor northerly azimuths
  const bool negate_latitude{!std::signbit(beta1.sin().v())};
  if (negate_latitude) {
    beta1 = -beta1;
    beta2 = -beta2;
  }

  // Ensure eastbound only
  const auto abs_lambda12{lambda12.abs()};

  // Estimate the azimuth at the start of the geodesic
  Angle<T> alpha0{
      (arc_length.v() >= ANTIPODAL_ARC_THRESHOLD)
          ? estimate_antipodal_initial_azimuth(beta1, beta2, abs_lambda12,
                                               ellipsoid)
      : estimate_azimuth
          ? estimate_initial_azimuth(beta1, beta2, abs_lambda12, ellipsoid)
          : great_circle::calculate_gc_azimuth(beta1, beta2, abs_lambda12)};

  // Use Newton's method to calculate the initial azimuth and aux length
  auto [alpha1, sigma12, iterations] =
      find_azimuth_arc_length_newtons_method<T>(
          beta1, beta2, abs_lambda12, alpha0, arc_length, tolerance, ellipsoid);

  // Calculate the correct azimuth for the start and finish points
  Angle<T> alpha2{calculate_end_azimuth(beta1, beta2, alpha1)};
  if (swap_latitudes) {
    std::swap(alpha1, alpha2);
  }

  if (swap_latitudes != negate_latitude) {
    alpha1 = alpha1.negate_cos();
    alpha2 = alpha2.negate_cos();
  }

  if (std::signbit(lambda12.sin().v())) {
    alpha1 = -alpha1;
    alpha2 = -alpha2;
  }

  Ensures((T() < sigma12.v()) && (sigma12.v() <= trig::PI<T>));

  return {alpha1, sigma12, alpha2, iterations};
}

/// Calculate the initial azimuth and great circle length between a pair
/// of points on the auxiliary sphere.
///
/// @pre beta_a.cos() >= 0
/// @pre beta_b.cos() >= 0
/// @pre -π <= delta_long <= π
/// @pre tolerance >= epsilon
///
/// @param beta1, beta2 the parametric latitudes of the start and finish
///     points on the auxiliary sphere.
/// @param delta_long the longitude difference.
/// @param tolerance the tolerance to perform the calculation to in Radians.
/// @param ellipsoid the `Ellipsoid`.
/// @return the azimuth at the start of the geodesic segment, the great circle
/// arc length on the auxiliary sphere, the azimuth at the end of
/// geodesic segment and the number of iterations required to calculate them.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto aux_sphere_azimuths_length(const Angle<T> beta1, const Angle<T> beta2,
                                const Angle<T> delta_long,
                                const Radians<T> tolerance,
                                const Ellipsoid<T> &ellipsoid)
    -> std::tuple<Angle<T>, Radians<T>, Angle<T>, unsigned> {
  const T MAX_EQUATORIAL_LENGTH{trig::PI<T> * ellipsoid.one_minus_f()};

  const Angle<T> gc_azimuth{
      great_circle::calculate_gc_azimuth(beta1, beta2, delta_long)};
  const auto gc_length{
      great_circle::calculate_gc_distance(beta1, beta2, delta_long)};

  // Determine whether on a meridian, i.e. a great circle which passes through
  // the North and South poles
  if (gc_azimuth.abs().sin().v() < std::numeric_limits<T>::epsilon()) {
    // gc_azimuth is 0° or 180°
    // Use opposite azimuth if points on opposite meridians
    const auto end_azimuth{(delta_long.cos().v() < T()) ? gc_azimuth.opposite()
                                                        : gc_azimuth};
    return {gc_azimuth, gc_length, end_azimuth, 0U};
  } else {
    // Determine whether on an equatorial path, i.e. the circle around the
    // equator.
    if ((gc_azimuth.cos().v() < std::numeric_limits<T>::epsilon()) &&
        (gc_length.v() < MAX_EQUATORIAL_LENGTH) &&
        (beta1.abs().sin().v() < std::numeric_limits<T>::epsilon()) &&
        (beta2.abs().sin().v() < std::numeric_limits<T>::epsilon())) {
      // Calculate the distance around the equator on the auxiliary sphere
      const Radians<T> equatorial_length{gc_length.v() *
                                         ellipsoid.recip_one_minus_f()};
      return {gc_azimuth, equatorial_length, gc_azimuth, 0U};
    } else {
      // Iterate to find the azimuth and length on the auxiliary sphere
      return find_azimuths_and_arc_length(beta1, beta2, delta_long, gc_length,
                                          tolerance, true, ellipsoid);
    }
  }
}

/// Calculate the `geodesic` azimuth and great circle length on the auxiliary
/// sphere between a pair of positions.
///
/// @pre a and b are valid LatLong's'
/// @pre tolerance >= epsilon
/// @post 0 <= arc_length <= PI
///
/// @param a, b the start and finish positions in geodetic coordinates.
/// @param tolerance the tolerance to perform the calculation to in Radians,
/// default great_circle::MIN_VALUE.
/// @param ellipsoid the `Ellipsoid`, default WGS 84.
/// @return the azimuth at the start of the geodesic segment, the great circle
/// arc length on the auxiliary sphere, the azimuth at the end of
/// geodesic segment and the number of iterations required to calculate them.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_azimuths_arc_length(
    const LatLong<T> &a, const LatLong<T> &b,
    const Radians<T> tolerance = Radians<T>(great_circle::MIN_VALUE<T>),
    const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
    -> std::tuple<Angle<T>, Radians<T>, Angle<T>, unsigned> {
  Expects(a.is_valid() && b.is_valid());

  // calculate the parametric latitudes on the auxiliary sphere
  const Angle<T> beta_a{
      ellipsoid.calculate_parametric_latitude(Angle<T>(a.lat()))};
  const Angle<T> beta_b{
      ellipsoid.calculate_parametric_latitude(Angle<T>(b.lat()))};

  // calculate the longitude difference
  const Angle<T> delta_long{Angle<T>(b.lon(), a.lon())};
  return aux_sphere_azimuths_length(beta_a, beta_b, delta_long, tolerance,
                                    ellipsoid);
}

/// Convert a great circle distance on the auxiliary sphere in radians to
/// metres on the ellipsoid.
///
/// @pre beta1.cos() >= 0
///
/// @param beta1 the start parametric Latitude on the auxiliary sphere.
/// @param alpha1 the azimuth at the start point.
/// @param arc_distance the great circle distance on the auxiliary sphere.
/// @param ellipsoid the `Ellipsoid`, default WGS 84.
///
/// @return the geodesic distance in Metres.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto convert_radians_to_metres(
    const Angle<T> &beta1, const Angle<T> &alpha1, Radians<T> arc_distance,
    const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84()) noexcept
    -> units::si::Metres<T> {

  Expects(beta1.cos().v() >= T());

  // Calculate the distance from the first equator crossing
  const auto sigma1{
      Angle<T>::from_y_x(beta1.sin().v(), beta1.cos().v() * alpha1.cos().v())};
  const auto sigma_sum{sigma1 + Angle<T>(arc_distance)};

  // Calculate the ellipsoid coefficients
  const trig::UnitNegRange<T> clairaut(alpha1.sin().v() * beta1.cos().v());
  const auto eps{ellipsoid.calculate_epsilon(clairaut)};
  const auto a1{ellipsoid::evaluate_A1(eps) + T(1)};
  const auto c1{ellipsoid::evaluate_coeffs_C1(eps)};
  const auto b11{ellipsoid::sin_cos_series(sigma1, c1)};
  const auto b12{ellipsoid::sin_cos_series(sigma_sum, c1)};

  return units::si::Metres<T>(ellipsoid.b().v() * a1 *
                              (arc_distance + b12 - b11).v());
}

} // namespace ellipsoid
} // namespace via
