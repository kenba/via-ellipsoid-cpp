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
/// @file ellipsoid.hpp
/// @brief Contains the via::ellipsoid software.
//////////////////////////////////////////////////////////////////////////////
/// @mainpage via-ellipsoid-cpp
///
/// A library for performing geometric calculations on the
/// [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System) ellipsoid.
///
/// This library uses the WGS84 primary parameters defined in Tab. 3-1 of the
/// [ICAO WGS 84 Implementation
/// Manual](https://www.icao.int/NACC/Documents/Meetings/2014/ECARAIM/REF08-Doc9674.pdf).
///
/// The shortest path between two points on the surface of an ellipsoid is a
/// [geodesic segment](https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid).
///
/// The equivalent of a straight line segment in planar geometry or a [great
/// circle arc](https://en.wikipedia.org/wiki/Great_circle) on the surface of a
/// sphere, see *Figure 1*.
///
/// <img
/// src="https://via-technology.aero/img/navigation/ellipsoid/sphere_mercator_long_geodesic.png"
/// width="600">
///
/// *Figure 1 A geodesic segment (orange) and great circle arc (blue)*
///
/// This library uses the correspondence between geodesic segments on an
/// ellipsoid and great-circle arcs on the auxiliary sphere, together with 3D
/// vectors to calculate:
///
/// - the length and azimuths of a geodesic segment between two positions;
/// - the along track and across track distances of a point relative
/// to a geodesic segment;
/// - and the intersection of two geodesic segments.
#include "ellipsoid/geodesic_intersection_functions.hpp"

namespace via {
namespace ellipsoid {

/// Calculate the position (Latitude and Longitude) where a pair of
/// `GeodesicSegment`s intersect, or None if the `GeodesicSegment`s do not
/// intersect.
/// @param g1, g2 the GeodesicSegments.
/// @param precision the precision in `Metres`.
///
/// @return the LatLong of the intersection point if the `GeodesicSegment`s
/// intersect, or std::nullopt if the points do not intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_intersection_point(const GeodesicSegment<T> &g1,
                                  const GeodesicSegment<T> &g2,
                                  units::si::Metres<T> precision)
    -> std::optional<LatLong<T>> {
  const Radians<T> precision_r{precision.v() / g1.ellipsoid().a().v()};
  const auto [distance1, distance2, _angle,
              _]{calculate_sphere_intersection_distances(g1, g2, precision_r)};
  if (vector::intersection::is_alongside(distance1, g1.arc_length(),
                                         precision_r) &&
      vector::intersection::is_alongside(distance2, g2.arc_length(),
                                         precision_r)) {
    // Ensure point is within g1
    const auto distance{distance1.clamp(g1.arc_length())};
    return g1.arc_lat_long(distance);
  }

  return std::nullopt;
}

} // namespace ellipsoid
} // namespace via
