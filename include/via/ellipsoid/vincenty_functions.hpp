#pragma once

//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2025 Ken Barker
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
/// @file vincenty_functions.hpp
/// @brief Contains the via::ellipsoid vincenty inverse template functions.
//////////////////////////////////////////////////////////////////////////////
#include "Ellipsoid.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>

namespace via {
namespace ellipsoid {
namespace vincenty {

/// Calculate the `geodesic` distance between a pair of positions.
///
/// Note: uses boost::geometry vincenty strategy.
/// @pre a and b are valid LatLong's'
///
/// @param a, b the start and finish positions in geodetic coordinates.
/// @param ellipsoid the `Ellipsoid`, default WGS 84.
/// @return the `geodesic` distance in `Metres`.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto inverse_distance(const LatLong<T> &a, const LatLong<T> &b,
                      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
    -> units::si::Metres<T> {
  using Wgs84Coords = boost::geometry::cs::geographic<boost::geometry::degree>;
  using GeographicPoint = boost::geometry::model::point<T, 2, Wgs84Coords>;

  // Note: the default boost geometry spheroid is WGS84
  // #include <boost/geometry/core/srs.hpp>
  using SpheroidType = boost::geometry::srs::spheroid<T>;

  // #include <boost/geometry/strategies/geographic/distance_vincenty.hpp>
  using VincentyStrategy =
      boost::geometry::strategy::distance::vincenty<SpheroidType>;

  SpheroidType spheriod(ellipsoid.a().v(), ellipsoid.b().v());
  VincentyStrategy vincenty(spheriod);

  // Note: order is longitude, latitude
  const GeographicPoint boost_a(a.lon().v(), a.lat().v());
  const GeographicPoint boost_b(b.lon().v(), b.lat().v());
  return units::si::Metres<T>(
      boost::geometry::distance(boost_a, boost_b, vincenty));
}

} // namespace vincenty
} // namespace ellipsoid
} // namespace via
