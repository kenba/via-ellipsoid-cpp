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
/// @file geodesic_intersection_functions.hpp
/// @brief Contains the via::ellipsoid geodesic intersection functions.
//////////////////////////////////////////////////////////////////////////////
#include "Geodesic.hpp"

namespace via {
namespace ellipsoid {
/// Clamp ref_length to zero if within precision or to arc_length if
/// within precision of arc_length
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto clamp_length(const T ref_length, const T arc_length,
                            const T precision) noexcept -> T {
  if (std::abs(ref_length) < precision)
    return T();
  else
    return (std::abs(ref_length - arc_length) < precision) ? arc_length
                                                           : ref_length;
}

/// Calculate the auxiliary Great Circle arc lengths to an intersection
/// point of two geodesics.
/// @param g1, g2 the Geodesics.
/// @param sq_precision the square of the Euclidean precision.
/// @param use_antipodal_intersection use the antipodal intersection point.
/// @param initial_distances the initial intersection distances in `Radians`.
/// @return the Great Circle lengths to the geodesic intersection point,
/// in radians.
template <typename T, unsigned MAX_ITERATIONS = 10u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_geodesic_intersection_distances(
    const Geodesic<T> &g1, const Geodesic<T> &g2, const T sq_precision,
    const bool use_antipodal_intersection,
    const std::tuple<Radians<T>, Radians<T>> initial_distances)
    -> std::tuple<Radians<T>, Radians<T>, unsigned> {
  auto [distance1, distance2]{initial_distances};

  auto iterations{1u};
  while (iterations < MAX_ITERATIONS) {
    const auto [pos1, pole1]{g1.aux_point_and_pole(distance1)};
    const auto [pos2, pole2]{g2.aux_point_and_pole(distance2)};

    const T sq_d{vector::sq_distance(pos1, pos2)};
    if (sq_d < sq_precision)
      break;

    ++iterations;

    // calculate the new intersection point
    const auto c{
        use_antipodal_intersection
            ? vector::intersection::calculate_intersection_point(pole2, pole1)
            : vector::intersection::calculate_intersection_point(pole1, pole2)};
    if (c.has_value()) {
      const auto [delta1, delta2]{
          vector::intersection::calculate_intersection_distances(
              pos1, pole1, pos2, pole2, c.value())};
      distance1 = distance1 + delta1;
      distance2 = distance2 + delta2;
    } else
      break;
  }

  return {distance1, distance2, iterations};
}

/// Calculate the distances along a pair of Geodesics (in Radians) to their
/// closest intersection or reference points.
/// @param g1, g2 the Geodesics.
/// @param precision the precision in `Radians`.
///
/// @return the distances along the Geodesics to the intersection point or to
/// their closest (reference) points if the Geodesics do not intersect and
/// the number of iterations required.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_aux_intersection_distances(const Geodesic<T> &g1,
                                          const Geodesic<T> &g2,
                                          Radians<T> precision)
    -> std::tuple<Radians<T>, Radians<T>, unsigned> {
  // The Geodesics MUST be on the same `Ellipsoid`

  // Convert precision in `Radians` to the square of Euclidean precision.
  const T e_precision{great_circle::gc2e_distance(precision)};
  const T sq_precision{e_precision * e_precision};

  // Get the start points and poles
  const auto [a1, pole1]{g1.aux_point_and_pole(Radians(T()))};
  const auto [a2, pole2]{g2.aux_point_and_pole(Radians(T()))};

  // if the start points are within precision of each other
  const T sq_d{vector::sq_distance(a1, a2)};
  if (sq_d < sq_precision)
    return {Radians(T()), Radians(T()), 0u};

  // Construct a geodesic between geodesic start points
  const auto [azi, aux_length]{aux_sphere_azimuth_length(
      g1.beta(), g2.beta(), g2.lon() - g1.lon(), g1.ellipsoid())};
  const Geodesic g3(g1.beta(), g1.lon(), azi, aux_length, g1.ellipsoid());

  // If the second geodesic start point lies on the first geodesic
  const auto [_a3, pole3]{g3.aux_point_and_pole(Radians(T()))};
  if (!vector::intersection::calculate_intersection_point(pole1, pole3)
           .has_value()) {
    const auto [atd, _xtd, iterations]{
        g1.calculate_aux_atd_and_xtd(g2.beta(), g2.lon(), precision)};
    // If the second geodesic end point lies on the first geodesic
    const auto [_a4, pole4]{g3.aux_point_and_pole(atd)};
    if (!vector::intersection::calculate_intersection_point(pole2, pole4)
             .has_value()) {
      // The geodesics are coincident
      const auto distances{
          vector::intersection::calculate_coincident_arc_distances(
              atd, pole1.dot(pole2) < T(), g1.aux_length(), g2.aux_length())};
      return {std::get<0>(distances), std::get<1>(distances), 0u};
    } else {
      // The geodesics intersect at the start of the second geodesic
      return {atd, Radians(0.0), iterations};
    }
  }

  // Determine whether the great circles on the auxiliary sphere are
  // coincident
  const auto c{
      vector::intersection::calculate_intersection_point(pole1, pole2)};
  if (c.has_value()) {
    const vector::Vector3<T> centroid{T(0.5) *
                                      (g1.mid_point() + g2.mid_point())};
    const bool use_antipodal_intersection =
        vector::intersection::use_antipodal_point(c.value(), centroid);
    const auto x{use_antipodal_intersection ? -c.value() : c.value()};
    const auto initial_distances{
        vector::intersection::calculate_intersection_distances(a1, pole1, a2,
                                                               pole2, x)};
    return calculate_geodesic_intersection_distances(
        g1, g2, sq_precision, use_antipodal_intersection, initial_distances);
  } else {
    const auto distances{
        vector::intersection::calculate_coincident_arc_distances(
            vector::calculate_great_circle_atd(a1, pole1, a2),
            pole1.dot(pole2) < T(), g1.aux_length(), g2.aux_length())};
    return {std::get<0>(distances), std::get<1>(distances), 0u};
  }
}

/// Calculate the distances along a pair of Geodesics (in Radians) to their
/// closest intersection or reference points.
/// @param g1, g2 the Geodesics.
/// @param precision the precision in `Metres`.
///
/// @return the distances along the Geodesics to the intersection point or to
/// their closest (reference) points if the Geodesics do not intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_intersection_distances(const Geodesic<T> &g1,
                                      const Geodesic<T> &g2,
                                      units::si::Metres<T> precision)
    -> std::tuple<Radians<T>, Radians<T>> {
  const Radians<T> precision_r{precision.v() / g1.ellipsoid().a().v()};
  const auto [distance1, distance2, _iterations]{
      calculate_aux_intersection_distances(g1, g2, precision_r)};
  return {distance1, distance2};
}

} // namespace ellipsoid
} // namespace via
