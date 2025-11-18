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
/// @file geodesic_intersection_functions.hpp
/// @brief Contains the via::ellipsoid geodesic segment intersection functions.
//////////////////////////////////////////////////////////////////////////////
#include "GeodesicSegment.hpp"
#include <cmath>
#include <via/sphere/great_circle.hpp>

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

/// Iterate the great circle arc distances to an intersection point of two
/// geodesic segments.
/// @param g1, g2 the GeodesicSegments.
/// @param sq_precision the square of the Euclidean precision.
/// @param use_antipodal_intersection use the antipodal intersection point.
/// @param initial_distances the initial intersection distances in `Radians`.
/// @return the Great Circle arc distances to the closest geodesic segment
/// intersection point, the angle at the intersection point and the
/// number of iterations required.
template <typename T, unsigned MAX_ITERATIONS = 10u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto iterate_geodesic_intersection_distances(
    const GeodesicSegment<T> &g1, const GeodesicSegment<T> &g2,
    const T sq_precision, const bool use_antipodal_intersection,
    const std::tuple<Radians<T>, Radians<T>> initial_distances)
    -> std::tuple<Radians<T>, Radians<T>, Angle<T>, unsigned> {
  auto [distance1, distance2]{initial_distances};

  auto iterations{1u};
  while (iterations < MAX_ITERATIONS) {
    const auto [pos1, pole1]{g1.arc_point_and_pole(distance1)};
    const auto [pos2, pole2]{g2.arc_point_and_pole(distance2)};

    const T sq_d{vector::sq_distance(pos1, pos2)};
    if (sq_d < sq_precision)
      break;

    ++iterations;

    // calculate the new intersection point
    const auto c{use_antipodal_intersection
                     ? vector::intersection::calculate_intersection(
                           pole2, pole1, vector::MIN_SQ_NORM<T>)
                     : vector::intersection::calculate_intersection(
                           pole1, pole2, vector::MIN_SQ_NORM<T>)};
    if (c.has_value()) {
      const auto [delta1, delta2]{
          vector::intersection::calculate_intersection_distances(
              pos1, pole1, pos2, pole2, c.value())};
      distance1 = distance1 + delta1;
      distance2 = distance2 + delta2;
    } else
      break;
  }

  // calculate the relative angle at the intersection point
  const auto angle{g2.arc_azimuth(Angle(distance2)) -
                   g1.arc_azimuth(Angle(distance1))};

  return {distance1, distance2, angle, iterations};
}

/// Calculate the distances along a pair of GeodesicSegments (in Radians)
/// to their closest intersection or reference points.
/// @pre g1.ellipsoid() == g2.ellipsoid()
/// @param g1, g2 the GeodesicSegments.
/// @param precision the precision in `Radians`.
///
/// @return the Great Circle arc distances to the closest geodesic segment
/// intersection point, the angle at the intersection point and the
/// number of iterations required.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_sphere_intersection_distances(const GeodesicSegment<T> &g1,
                                             const GeodesicSegment<T> &g2,
                                             Radians<T> precision)
    -> std::tuple<Radians<T>, Radians<T>, Angle<T>, unsigned> {
  // The GeodesicSegments MUST be on the same `Ellipsoid`
  Expects(g1.ellipsoid() == g2.ellipsoid());

  // Convert precision in `Radians` to the square of Euclidean precision.
  const T e_precision{great_circle::gc2e_distance(precision)};
  const T sq_precision{e_precision * e_precision};

  // Determine whether the geodesics are reciprocal
  const Angle<T> delta_azimuth1_2{g2.azi() - g1.azi()};
  const bool reciprocal{std::signbit(delta_azimuth1_2.cos().v())};

  // if the start points are within precision of each other
  const T sq_d{vector::sq_distance(g1.a(), g2.a())};
  if (sq_d < sq_precision) {
    return {Radians(T()), Radians(T()), delta_azimuth1_2, 0u};
  }

  // Construct a geodesic between geodesic start points
  const auto [g3_azi, g3_arc_length, g3_end_azi, _]{aux_sphere_azimuths_length(
      g1.beta(), g2.beta(), g2.lon() - g1.lon(),
      Radians<T>(great_circle::MIN_VALUE<T>), g1.ellipsoid())};

  const Radians<T> atd{reciprocal ? -g3_arc_length : g3_arc_length};

  // Determine whether the geodesics are coincident
  const Angle<T> delta_azimuth1_3{g3_azi - g1.azi()};
  const Angle<T> delta_azimuth2_3{g3_end_azi - g2.azi()};
  if ((delta_azimuth1_3.sin().abs().v() < vector::MIN_SIN_ANGLE<T>) &&
      (delta_azimuth2_3.sin().abs().v() < vector::MIN_SIN_ANGLE<T>)) {
    // The geodesics are coincident
    const auto [distance1, distance2]{
        vector::intersection::calculate_coincident_arc_distances(
            atd, reciprocal, g1.arc_length(), g2.arc_length())};
    const auto angle{reciprocal ? Angle<T>().opposite() : Angle<T>()};
    return {distance1, distance2, angle, 0u};
  }

  // Calculate the intersection of the poles at the mid points of the unit
  // sphere great circle arcs
  const Radians<T> half_arc_length1{g1.arc_length().half()};
  const Radians<T> half_arc_length2{g2.arc_length().half()};
  const auto [a1mid, pole1mid]{g1.arc_point_and_pole(half_arc_length1)};
  const auto [a2mid, pole2mid]{g2.arc_point_and_pole(half_arc_length2)};
  const auto c{vector::intersection::calculate_intersection(
      pole1mid, pole2mid, vector::MIN_SQ_NORM<T>)};

  // Determine whether the great circles on the auxiliary sphere are
  // coincident
  if (c.has_value()) {
    // intersesction found, chose the closest intersesction point to the
    // centroid of the unit sphere great circle arc mid points
    const vector::Vector3<T> centroid{T(0.5) * (a1mid + a2mid)};
    const bool use_antipodal_intersection =
        vector::intersection::use_antipodal_point(c.value(), centroid);
    const auto x{use_antipodal_intersection ? -c.value() : c.value()};

    auto [d1, d2]{vector::intersection::calculate_intersection_distances(
        a1mid, pole1mid, a2mid, pole2mid, x)};
    d1 += half_arc_length1;
    d2 += half_arc_length2;

    return iterate_geodesic_intersection_distances(
        g1, g2, sq_precision, use_antipodal_intersection, {d1, d2});
  } else {
    // This code should never be executed.
    // The check for coincident geodesics should cover coincident great circles.
    const auto [distance1, distance2]{
        vector::intersection::calculate_coincident_arc_distances(
            atd, reciprocal, g1.arc_length(), g2.arc_length())};
    const auto angle{reciprocal ? Angle<T>().opposite() : Angle<T>()};
    return {distance1, distance2, angle, 0u};
  }
}

/// Calculate the distances along a pair of GeodesicSegments (in Radians) to
/// their closest intersection or reference points.
/// @param g1, g2 the GeodesicSegments.
/// @param precision the precision in `Metres`.
///
/// @return the distances along the GeodesicSegments to the intersection point
/// or to their closest (reference) points if the GeodesicSegments do not
/// intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_intersection_distances(const GeodesicSegment<T> &g1,
                                      const GeodesicSegment<T> &g2,
                                      units::si::Metres<T> precision)
    -> std::tuple<Radians<T>, Radians<T>> {
  const Radians<T> precision_r{precision.v() / g1.ellipsoid().a().v()};
  const auto [distance1, distance2, _angle, _iterations]{
      calculate_sphere_intersection_distances(g1, g2, precision_r)};
  return {distance1, distance2};
}

} // namespace ellipsoid
} // namespace via
