#pragma once

//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2019-2026 Ken Barker
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
#include <via/angle.hpp>
#include <via/sphere/great_circle.hpp>
#include <via/sphere/vector.hpp>

namespace via {
namespace ellipsoid {
namespace intersection {

/// Determine whether two `GeodesicSegment`s are coincident,
/// i.e. they lie on the same geodesic path.
///
/// @param g_0, g_1 the `GeodesicSegment`s.
/// @param min_sin_angle the sine of the minimum angle between coincident
/// `GeodesicSegment`s.
///
/// @return true if the segments are on the same geodesic path, false otherwise.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto geodesics_are_coincident(const GeodesicSegment<T> &g_0,
                              const GeodesicSegment<T> &g_1, T min_sin_angle)
    -> bool {
  // construct a geodesic between geodesic start points
  const auto [g_2_azi, _g_2_arc_length, g_2_end_azi, _]{
      aux_sphere_azimuths_length(g_0.beta(), g_1.beta(), g_1.lon() - g_0.lon(),
                                 Radians<T>(great_circle::MIN_VALUE<T>),
                                 g_0.ellipsoid())};

  // the segments are coincident if their angle differences to the
  // path between start positions are both within min_sin_angle
  const Angle<T> delta_azimuth0_2{g_2_azi - g_0.azi()};
  const Angle<T> delta_azimuth1_2{g_2_end_azi - g_1.azi()};

  return (delta_azimuth0_2.sin().abs().v() < min_sin_angle) &&
         (delta_azimuth1_2.sin().abs().v() < min_sin_angle);
}

/// Find the closest intersection distances of two `GeodesicSegment`s.
///
/// @param g_0, g_1 the `GeodesicSegment`s.
/// @param use_antipodal_intersection use the antipodal intersection point
/// @param distance_0, distance_1 the initial arc intersection distances from
/// the segment start.
/// @param precision the precision in `Radians`
///
/// @return the arc distances from the segment start to the closest intersection
/// point, the relative angle at the intersection point, and the number of
/// iterations required.
template <typename T, unsigned MAX_ITERATIONS = 10u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto find_geodesic_intersection_distances(const GeodesicSegment<T> &g_0,
                                          const GeodesicSegment<T> &g_1,
                                          const bool use_antipodal_intersection,
                                          Radians<T> distance_0,
                                          Radians<T> distance_1,
                                          const Radians<T> precision)
    -> std::tuple<Radians<T>, Radians<T>, unsigned> {
  const auto sq_precision{std::pow(great_circle::gc2e_distance(precision), 2)};

  auto iterations{1u};
  while (iterations < MAX_ITERATIONS) {
    const auto [point_0, pole_0]{g_0.arc_point_and_pole(distance_0)};
    const auto [point_1, pole_1]{g_1.arc_point_and_pole(distance_1)};

    const T sq_d{vector::sq_distance(point_0, point_1)};
    if (sq_d < sq_precision)
      break;

    ++iterations;

    // calculate the new intersection point
    const auto x{use_antipodal_intersection
                     ? vector::intersection::calculate_intersection(
                           pole_1, pole_0, vector::MIN_SQ_NORM<T>)
                     : vector::intersection::calculate_intersection(
                           pole_0, pole_1, vector::MIN_SQ_NORM<T>)};
    if (x.has_value()) {
      distance_0 += vector::calculate_great_circle_atd(point_0, pole_0, *x);
      distance_1 += vector::calculate_great_circle_atd(point_1, pole_1, *x);
    } else
      break;
  }

  return {distance_0, distance_1, iterations};
}

/// Find the closest intersection distances of two `GeodesicSegment`s and the
/// relative angle at the reference point.
///
/// The reference point is the closest intersection point or the centroid
/// if the `GeodesicSegment`s are coincident.
///
/// @param g_0, g_1 the `GeodesicSegment`s.
/// @param precision the precision in `Radians`.
///
/// @return the arc distances from the segment mid points to the reference
/// point, the relative angle at the reference point, and the number of
/// iterations required.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_arc_reference_distances_and_angle(const GeodesicSegment<T> &g_0,
                                                 const GeodesicSegment<T> &g_1,
                                                 Radians<T> precision)
    -> std::tuple<Radians<T>, Radians<T>, Angle<T>, unsigned> {
  // The GeodesicSegments MUST be on the same `Ellipsoid`
  Expects(g_0.ellipsoid() == g_1.ellipsoid());

  if (g_0 == g_1) {
    return {Radians(T()), Radians(T()), Angle<T>(), T()};
  }

  const Radians<T> half_length_0{g_0.arc_length().half()};
  const Radians<T> half_length_1{g_1.arc_length().half()};
  const auto [mid_point_0, pole_0]{g_0.arc_point_and_pole(half_length_0)};
  const auto [mid_point_1, pole_1]{g_1.arc_point_and_pole(half_length_1)};
  const vector::Vector3<T> centroid{(mid_point_0 + mid_point_1)};
  const auto intersection{vector::intersection::calculate_intersection(
      pole_0, pole_1, vector::MIN_SQ_NORM<T>)};
  if (intersection.has_value() &&
      !geodesics_are_coincident(g_0, g_1, vector::MIN_SIN_ANGLE<T>)) {
    // great circles or geodesics interact

    // find the closest intersection
    const vector::Vector3<T> c{intersection.value()};
    const bool use_antipodal_intersection =
        vector::intersection::use_antipodal_point(c, centroid);
    const vector::Vector3<T> x{use_antipodal_intersection ? -c : c};

    // calculate distances to the closest geodesic intersection from arc start
    const auto [distance_0, distance_1,
                iterations]{find_geodesic_intersection_distances(
        g_0, g_1, use_antipodal_intersection,
        vector::calculate_great_circle_atd(g_0.a(), pole_0, x),
        vector::calculate_great_circle_atd(g_1.a(), pole_1, x), precision)};

    const Angle<T> angle{g_1.arc_azimuth(Angle(distance_1)) -
                         g_0.arc_azimuth(Angle(distance_0))};

    return {distance_0 - half_length_0, distance_1 - half_length_1, angle.abs(),
            iterations};
  } else {
    // great circles or geodesics are coincident

    const vector::Vector3<T> c{
        vector::normalise_centroid(centroid, mid_point_0, pole_0)};
    const Radians<T> distance_0{
        vector::calculate_great_circle_atd(mid_point_0, pole_0, c)};
    const Radians<T> distance_1{
        vector::calculate_great_circle_atd(mid_point_1, pole_1, c)};

    const Angle<T> angle{
        std::signbit(pole_0.dot(pole_1)) ? Angle<T>().opposite() : Angle<T>()};

    return {distance_0, distance_1, angle, T()};
  }
}

} // namespace intersection
} // namespace ellipsoid
} // namespace via
