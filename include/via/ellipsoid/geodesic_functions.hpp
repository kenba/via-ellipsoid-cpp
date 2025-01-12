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
/// @file geodesic_functions.hpp
/// @brief Contains the via::ellipsoid geodesic template functions.
//////////////////////////////////////////////////////////////////////////////
#include "Ellipsoid.hpp"
#include <via/angle/trig.hpp>
#include <via/sphere.hpp>
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
  if ((q <= T()) && (r <= T()))
    return T();

  const T S{p * q / T(4)};
  const T r2{r * r};
  const T r3{r * r2};
  T u{r};

  // The discriminant of the quadratic equation for T3.
  // This is zero on the evolute curve p^(1/3)+q^(1/3) = 1
  const T discriminant{S * (S + 2 * r3)};
  if (T() <= discriminant) {
    T t3{S + r3};
    // Pick the sign on the sqrt to maximize abs(T3), to minimise loss
    // of precision due to cancellation.
    t3 += std::copysign(std::sqrt(discriminant), t3);
    const T t{std::cbrt(t3)};
    u += t + (t != T() ? r2 / t : T());
  } else {
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

/// Calculate: m12b = (reduced length)/_b
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_reduced_length(const T eps, const Radians<T> sigma12,
                                        const Angle<T> sigma1, const T dn1,
                                        const Angle<T> sigma2, const T dn2)
    -> T {
  // TODO investigate why this precondition fails
  // Expects(T() <= sigma12.v());

  const T A1{evaluate_A1<T>(eps)};
  const T A2{evaluate_A2<T>(eps)};
  const T m0x{A1 - A2};

  const T a1p1 = T(1) + A1;
  const T a2p1 = T(1) + A2;

  const auto Ca{evaluate_coeffs_C1<T>(eps)};
  auto Cb{evaluate_coeffs_C2<T>(eps)};

  // Assume here that Ca.size() >= Cb.size()
  // int size(Cb.size() - 1);
  for (auto i(1u); i < Cb.size(); ++i)
    Cb[i] = a1p1 * Ca[i] - a2p1 * Cb[i];

  const T J12{m0x * sigma12.v() + (sin_cos_series(sigma2, Cb).v() -
                                   sin_cos_series(sigma1, Cb).v())};
  return dn2 * (sigma1.cos().v() * sigma2.sin().v()) -
         dn1 * (sigma1.sin().v() * sigma2.cos().v()) -
         sigma1.cos().v() * sigma2.cos().v() * J12;
}

/// Estimate the initial azimuth on the auxiliary sphere for a nearly antipodal
/// arc. It calculates and solves the astroid problem.
/// @pre 0 <= lambda12.sin()
/// @param beta1, beta2 the parametric latitudes of the start and
/// finish points on the auxiliary sphere.
/// @param lambda12 longitude difference between start and finish points.
/// @param ellipsoid the `Ellipsoid`.
/// @return the estimate of the initial azimuth on the auxiliary sphere.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto estimate_antipodal_initial_azimuth(const Angle<T> beta1,
                                                  const Angle<T> beta2,
                                                  const Angle<T> lambda12,
                                                  const Ellipsoid<T> &ellipsoid)
    -> Angle<T> {
#ifdef _MSC_VER
  const T X_THRESHOLD{1000 * std::sqrt(std::numeric_limits<T>::epsilon())};
#else
  constexpr T X_THRESHOLD{1000 * std::sqrt(std::numeric_limits<T>::epsilon())};
#endif
  constexpr T Y_TOLERANCE{200 * std::numeric_limits<T>::epsilon()};

  Expects(T() <= lambda12.sin().v());

  // Calculate the integration parameter for geodesic
  const auto clairaut{beta1.cos()}; // Note: assumes sin_alpha_1 = 1
  const T eps{ellipsoid.calculate_epsilon(clairaut)};
  const T a3f{evaluate_poynomial(eps, ellipsoid.a3())};

  const T lamscale{T(ellipsoid.f()) * beta1.cos().v() * a3f * trig::PI<T>};
  const T betscale{lamscale * beta1.cos().v()};

  // Solve astroid problem
  const T x{lambda12.opposite().to_radians().v() / lamscale};
  const T y{sine_sum(beta1, beta2).v() / betscale};

  // Test x and y params
  if ((y > -Y_TOLERANCE) && (x > T(-1) - X_THRESHOLD)) {
    const trig::UnitNegRange<T> sin_alpha{std::min(-x, T(1))};
    return Angle<T>(sin_alpha, trig::cosine_from_sine(sin_alpha, T(-1)));
  } else {
    const T k{calculate_astroid(x, y)};
    const Radians<T> omg12a{lamscale * (-x * k / (1 + k))};

    Angle<T> omega12(omg12a);
    return great_circle::calculate_gc_azimuth(beta1, beta2,
                                              omega12.negate_cos());
  }
}

/// Calculate the cosine of the longitude difference from the equator crossing.
/// @param beta the `parametric` latitude
/// @param cos_azimuth the cosine of the azimuth at the `parametric` latitude
///
/// @return the cosine of the longitude difference, zero if the `parametric`
/// latitude is very close to the equator.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_cos_omega(const Angle<T> beta,
                                   const trig::UnitNegRange<T> cos_azimuth)
    -> T {
  return (beta.sin().abs().v() < std::numeric_limits<T>::epsilon())
             ? T(1)
             : cos_azimuth.v() * beta.cos().v();
}

/// Calculate the azimuth on the auxiliary sphere at latitude beta2.
/// @pre 0 <= alpha1.sin()
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto calculate_end_azimuth(const Angle<T> beta1, const Angle<T> beta2,
                                     const Angle<T> alpha1) -> Angle<T> {
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
    const T temp4{(T() < temp3) ? std::sqrt(temp3) / beta2.cos().v() : T()};
    cos_alpha2 = trig::UnitNegRange<T>::clamp(temp4);
  }

  return Angle<T>{sin_alpha2, cos_alpha2};
}

/// Calculate the longitude difference between the auxiliary sphere and
/// ellipsoid.
/// @pre 0 < sigma12
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
  return Radians<T>(a3c * (sigma12.v() + (b32.v() - b31.v())));
}

/// Find the aziumth and great circle length on the auxiliary sphere.
/// It uses Newton's method to solve:
///   f(alp1) = lambda12(alp1) - lam12 = 0
/// @tparam MAX_ITERS the maximum numer of iterations to attempt.
/// @param beta_a, beta_b the geodetic latitudes of the start and finish points.
/// @param lambda12 the geodesic longitude difference in radians.
/// @param gc_length the auxiliary sphere great circle length in radians
/// @param ellipsoid the `Ellipsoid`.
/// @param tolerance the tolerance to perform the calculation to in Radians.
/// @return the azimuth and great circle length on the auxiliary sphere at the
/// start of the geodesic and the number of iterations required to calculate
/// them.
template <typename T, int MAX_ITERS = 20>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto find_azimuth_and_aux_length(const Angle<T> beta_a, const Angle<T> beta_b,
                                 const Angle<T> lambda12,
                                 const Radians<T> gc_length,
                                 const Ellipsoid<T> &ellipsoid,
                                 const Radians<T> tolerance)
    -> std::tuple<Angle<T>, Radians<T>, unsigned> {
  const T ANTIPODAL_ARC_THRESHOLD{trig::PI<T> * ellipsoid.one_minus_f()};
  Expects(T() < gc_length.v());

  // Start at the latitude furthest from the Equator
  const bool swap_latitudes{beta_a.sin().abs() < beta_b.sin().abs()};
  Angle<T> beta1{swap_latitudes ? beta_b : beta_a};
  Angle<T> beta2{swap_latitudes ? beta_a : beta_b};

  // Start South of the Equator
  const bool negate_latitude{T() < beta1.sin().v()};
  if (negate_latitude) {
    beta1 = -beta1;
    beta2 = -beta2;
  }

  const T dn1{std::sqrt(T(1) + T(ellipsoid.ep_2()) * beta1.sin().v() *
                                   beta1.sin().v())};
  const T dn2{std::sqrt(T(1) + T(ellipsoid.ep_2()) * beta2.sin().v() *
                                   beta2.sin().v())};

  // Use positive longitude difference.
  const Angle<T> abs_lambda12{lambda12.abs()};

  // Estimate the azimuth at the start of the geodesic
  Angle<T> alpha1(Radians<T>(0));
  if (ANTIPODAL_ARC_THRESHOLD < gc_length.v()) {
    alpha1 = estimate_antipodal_initial_azimuth(beta1, beta2, abs_lambda12,
                                                ellipsoid);
  } else {
    alpha1 = great_circle::calculate_gc_azimuth(beta1, beta2, abs_lambda12);
  }
  Angle<T> alpha2{alpha1};
  auto sigma12{gc_length};

#ifdef OUTPUT_GEOD_ITERATOR_STEPS
  auto prev_v{gc_length.v()};
#endif
  auto iterations{1U};
  for (; iterations <= MAX_ITERS; ++iterations) {
    // Calculate Clairaut's constant
    const trig::UnitNegRange<T> clairaut(alpha1.sin().v() * beta1.cos().v());
    const T eps{ellipsoid.calculate_epsilon(clairaut)};

    // Calculate first longitude (omega1) and distance (sigma1) from the
    // Northbound equator crossing
    const trig::UnitNegRange<T> sin_omega1{clairaut.v() * beta1.sin().v()};
    const auto cos_omega1{calculate_cos_omega(beta1, alpha1.cos())};
    const auto omega1{Angle<T>::from_y_x(sin_omega1.v(), cos_omega1)};
    const auto sigma1{Angle<T>::from_y_x(beta1.sin().v(), cos_omega1)};

    // Calculate azimuth at the end point
    alpha2 = calculate_end_azimuth(beta1, beta2, alpha1);

    // Calculate second longitude (omega2) and distance (sigma2) from the
    // Northbound equator crossing
    const trig::UnitNegRange<T> sin_omega2{clairaut.v() * beta2.sin().v()};
    const auto cos_omega2{calculate_cos_omega(beta2, alpha2.cos())};
    const auto omega2{Angle<T>::from_y_x(sin_omega2.v(), cos_omega2)};
    const auto sigma2{Angle<T>::from_y_x(beta2.sin().v(), cos_omega2)};

    // Calculate Longitude difference on the auxiliary sphere
    auto omega12 = omega2 - omega1;
    // clamp to range 0 to Pi
    if (omega12.sin().v() < T())
      omega12 = (omega12.cos().v() < T()) ? Angle<T>(trig::UnitNegRange(T()),
                                                     trig::UnitNegRange(T(-1)))
                                          : Angle<T>();

    // Calculate great circle length on the auxiliary sphere
    Angle<T> sc_sigma12{sigma2 - sigma1};
    if (sc_sigma12.sin().v() < T()) // clamp to range 0 to Pi
      sc_sigma12 =
          (sc_sigma12.cos().v() < T())
              ? Angle<T>(trig::UnitNegRange(T()), trig::UnitNegRange(T(-1)))
              : Angle<T>();
    sigma12 = sc_sigma12.abs().to_radians();

    // Calculate difference between geodesic and great circle longitudes
    const auto eta{(omega12 - abs_lambda12).to_radians()};
    const auto domg12{
        delta_omega12(clairaut, eps, sigma12, sigma1, sigma2, ellipsoid)};
    const T v{eta.v() - domg12.v()};

    // Test within tolerance
    if (std::abs(v) <= tolerance.v())
      break;

    // Calculate the denominator for Newton's method
    const T dv{(alpha2.cos().abs().v() <= std::numeric_limits<T>::epsilon())
                   ? -2 * ellipsoid.one_minus_f() * dn1 / beta1.sin().v()
                   : ellipsoid.one_minus_f() *
                         calculate_reduced_length(eps, sigma12, sigma1, dn1,
                                                  sigma2, dn2) /
                         (alpha2.cos().v() * beta2.cos().v())};

    // Calculate the change in initial azimuth and test within tolerance
    const T dalpha1{std::clamp<T>(-v / dv, -1, 1)};
    if (std::abs(dalpha1) <= tolerance.v())
      break;

#ifdef OUTPUT_GEOD_ITERATOR_STEPS
    const bool converging{std::abs(v) < prev_v};
    prev_v = std::abs(v);

    std::cout << std::boolalpha;
    std::cout << i << ',' << converging << ',' << dalpha1 << ',' << v << ','
              << dv << std::endl;
#endif

    // Adjust the azimuth by dalpha1 for the next iteration
    alpha1 = alpha1 + Angle<T>(Radians<T>{dalpha1});
  }

  if (swap_latitudes)
    alpha1 = alpha2;

  if (swap_latitudes != negate_latitude)
    alpha1 = alpha1.negate_cos();

  const bool lambda12_negative{lambda12.sin().v() < T()};
  if (lambda12_negative)
    alpha1 = -alpha1;

  Ensures(T() < sigma12.v());

  return {alpha1, sigma12, iterations};
}

/// Calculate the initial azimuth and great circle length between a pair
/// of points on the auxiliary sphere.
/// @param beta1, beta2 the parametric latitudes of the start and finish
///     points on the auxiliary sphere.
/// @param delta_long the longitude difference.
/// @param ellipsoid the `Ellipsoid`.
/// @param tolerance the tolerance to perform the calculation to in Radians.
/// @return the azimuth and great circle length on the auxiliary sphere at the
/// start of the geodesic and the number of iterations required to calculate
/// them.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto aux_sphere_azimuth_length(const Angle<T> beta1, const Angle<T> beta2,
                               const Angle<T> delta_long,
                               const Ellipsoid<T> &ellipsoid,
                               const Radians<T> tolerance)
    -> std::tuple<Angle<T>, Radians<T>, unsigned> {
  const Angle<T> gc_azimuth{
      great_circle::calculate_gc_azimuth(beta1, beta2, delta_long)};
  const auto gc_length{
      great_circle::calculate_gc_distance(beta1, beta2, delta_long)};

  // Determine whether on a meridian, i.e. a great circle which passes through
  // the North and South poles
  if (gc_azimuth.abs().sin().v() < great_circle::MIN_VALUE<T>) {
    // gc_azimuth is 0° or 180°
    return {gc_azimuth, gc_length, 0U};
  } else {
    // Determine whether on an equatorial path, i.e. the circle around the
    // equator.
    if ((gc_azimuth.cos().v() < great_circle::MIN_VALUE<T>) &&
        (beta1.abs().sin().v() < std::numeric_limits<T>::epsilon()) &&
        (beta2.abs().sin().v() < std::numeric_limits<T>::epsilon())) {
      // Calculate the distance around the equator on the auxiliary sphere
      const Radians<T> equatorial_length{gc_length.v() *
                                         ellipsoid.recip_one_minus_f()};
      return {gc_azimuth, equatorial_length, 0U};
    } else {
      // Iterate to find the azimuth and length on the auxiliary sphere
      return find_azimuth_and_aux_length(beta1, beta2, delta_long, gc_length,
                                         ellipsoid, tolerance);
    }
  }
}

/// Calculate the `geodesic` azimuth and great circle length on the auxiliary
/// sphere between a pair of positions.
/// @pre a and b are valid LatLong's'
/// @post 0 <= aux_length <= PI
/// @param a, b the start and finish positions in geodetic coordinates.
/// @param ellipsoid the `Ellipsoid`.
/// @param tolerance the tolerance to perform the calculation to in Radians.
/// @return the azimuth and great circle length on the auxiliary sphere at the
/// start of the geodesic and the number of iterations required to calculate
/// them.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_azimuth_aux_length(
    const LatLong<T> &a, const LatLong<T> &b, const Ellipsoid<T> &ellipsoid,
    const Radians<T> tolerance = Radians<T>(great_circle::MIN_VALUE<T>))
    -> std::tuple<Angle<T>, Radians<T>, unsigned> {
  Expects(a.is_valid() && b.is_valid());

  // calculate the parametric latitudes on the auxiliary sphere
  const Angle<T> beta_a{
      ellipsoid.calculate_parametric_latitude(Angle<T>(a.lat()))};
  const Angle<T> beta_b{
      ellipsoid.calculate_parametric_latitude(Angle<T>(b.lat()))};

  // calculate the longitude difference
  const Angle<T> delta_long{Angle<T>(b.lon(), a.lon())};
  return aux_sphere_azimuth_length(beta_a, beta_b, delta_long, ellipsoid,
                                   tolerance);
}

/// Convert a great circle distance on the auxiliary sphere in radians to
/// metres on the ellipsoid.
/// @param beta1 the start parametric Latitude on the auxiliary sphere.
/// @param alpha1 the azimuth at the start point.
/// @param gc_distance the great circle distance on the auxiliary sphere in
/// radians.
/// @param ellipsoid the `Ellipsoid`.
///
/// @return the geodesic distance in Metres.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto convert_radians_to_metres(const Angle<T> &beta1,
                                         const Angle<T> &alpha1,
                                         Radians<T> gc_distance,
                                         const Ellipsoid<T> &ellipsoid) noexcept
    -> units::si::Metres<T> {
  // Calculate the distance from the first equator crossing
  const auto cos_omega1{calculate_cos_omega(beta1, alpha1.cos())};
  const auto sigma1{Angle<T>::from_y_x(beta1.sin().v(), cos_omega1)};
  const auto sigma_sum{sigma1 + Angle<T>(gc_distance)};

  // Calculate the ellipsoid coefficients
  const trig::UnitNegRange<T> clairaut(alpha1.sin().v() * beta1.cos().v());
  const auto eps{ellipsoid.calculate_epsilon(clairaut)};
  const auto a1{ellipsoid::evaluate_A1(eps) + T(1)};
  const auto c1{ellipsoid::evaluate_coeffs_C1(eps)};
  const auto b11{ellipsoid::sin_cos_series(sigma1, c1)};
  const auto b12{ellipsoid::sin_cos_series(sigma_sum, c1)};

  return units::si::Metres<T>(ellipsoid.b().v() * a1 *
                              (gc_distance + b12 - b11).v());
}

} // namespace ellipsoid
} // namespace via
