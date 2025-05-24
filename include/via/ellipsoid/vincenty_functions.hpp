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
#include <boost/geometry/formulas/vincenty_inverse.hpp>
#include <boost/geometry/srs/spheroid.hpp>

namespace via {
namespace ellipsoid {
namespace vincenty {

/// Calculate the `geodesic` azimuths and distance between a pair of positions.
///
/// Note: uses boost::geometry vincenty_inverse formula.
/// @pre a and b are valid LatLong's'
///
/// @param a, b the start and finish positions in geodetic coordinates.
/// @param ellipsoid the `Ellipsoid`, default WGS 84.
/// @return the azimuth at the start of the geodesic segment,
/// the `geodesic` distance in `Metres` and the azimuth at the end of
/// the geodesic segment.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto inverse_azimuths_and_distance(
    const LatLong<T> &a, const LatLong<T> &b,
    const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
    -> std::tuple<Radians<T>, units::si::Metres<T>, Radians<T>> {
  using namespace boost::geometry;

  const auto result{formula::vincenty_inverse<T, true, true, true>().apply(
      trig::deg2rad(a.lon().v()), trig::deg2rad(a.lat().v()),
      trig::deg2rad(b.lon().v()), trig::deg2rad(b.lat().v()),
      srs::spheroid<T>(ellipsoid.a().v(), ellipsoid.b().v()))};

  return {Radians<T>(result.azimuth), units::si::Metres<T>(result.distance),
          Radians<T>(result.reverse_azimuth)};
}

} // namespace vincenty
} // namespace ellipsoid
} // namespace via
