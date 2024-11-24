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
/// @file ellipsoid.hpp
/// @brief Contains the via::ellipsoid software.
//////////////////////////////////////////////////////////////////////////////
#include "ellipsoid/geodesic_intersection_functions.hpp"

namespace via {
namespace ellipsoid {

/// Calculate the position (Latitude and Longitude) where a pair of
/// `Geodesic`s intersect, or None if the `Geodesic`s do not intersect.
/// @param g1, g2 the Geodesics.
/// @param precision the precision in `Metres`.
///
/// @return the LatLong of the intersection point if the points intersect,
/// or std::nullopt if the points do not intersect.
template <typename T>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
auto calculate_intersection_point(const Geodesic<T> &g1, const Geodesic<T> &g2,
                                  units::si::Metres<T> precision)
    -> std::optional<LatLong<T>> {
  const auto [distance1,
              distance2]{calculate_intersection_distances(g1, g2, precision)};
  if (vector::intersection::is_within(distance1.v(), g1.aux_length().v()) &&
      vector::intersection::is_within(distance2.v(), g2.aux_length().v()))
    return g1.aux_lat_long(distance1);
  else
    return std::nullopt;
}

} // namespace ellipsoid
} // namespace via
