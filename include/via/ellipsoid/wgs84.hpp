#pragma once

//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2016-2024 Ken Barker
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
/// @file wgs84.hpp
/// @brief Contains the primary parameters of the WGS84 ellipsoid.
/// From Eurocontrol [WGS 84 Implementation Manual
/// Version 2.4](https://www.icao.int/safety/pbn/Documentation/EUROCONTRL/Eurocontrol%20WGS%2084%20Implementation%20Manual.pdf)
/// Chapter 3, page 14.
#include "Metres.hpp"

namespace via {
namespace ellipsoid {
namespace wgs84 {
/// The WGS 84 Semimajor axis measured in metres.
/// This is the radius at the equator.
template <typename T>
  requires std::floating_point<T>
constexpr Metres<T> A{6'378'137};

/// The WGS 84 flattening, a ratio.
/// This is the flattening of the ellipsoid at the poles.
template <typename T>
  requires std::floating_point<T>
constexpr T F{static_cast<T>(1.0L / 298.257'223'563L)};
} // namespace wgs84
} // namespace ellipsoid
} // namespace via
