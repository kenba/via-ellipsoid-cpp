#!/usr/bin/env python

# Copyright (c) 2019-2025 Ken Barker
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#  @file test_Geodesic
#  @brief Contains unit tests for the via ellipsoid Geodesic class.

import pytest
from numpy.testing import assert_almost_equal
from via_angle import Degrees
from via_sphere import LatLong
from via_units import Metres
from via_ellipsoid import Ellipsoid, Geodesic, calculate_intersection_point

def test_intersection_point_distance():
    # Karney's example:
    # Istanbul, Washington, Reyjavik and Accra
    # from:
    # <https://sourceforge.net/p/geographiclib/discussion/1026621/thread/21aaff9f/#fe0a>
    istanbul = LatLong(Degrees(42.0), Degrees(29.0))
    washington = LatLong(Degrees(39.0), Degrees(-77.0))
    reyjavik = LatLong(Degrees(64.0), Degrees(-22.0))
    accra = LatLong(Degrees(6.0), Degrees(0.0))

    g1 = Geodesic(istanbul, washington)
    g2 = Geodesic(reyjavik, accra)

    intersection_point_1 = calculate_intersection_point(g1, g2, Metres(1e-3))
    if intersection_point_1:
        assert_almost_equal(54.7170296089477, intersection_point_1.lat().v())
        assert_almost_equal(-14.56385574430775, intersection_point_1.lon().v())

    # Swap Geodesics
    intersection_point_2 = calculate_intersection_point(g2, g1, Metres(1e-3))
    if intersection_point_2:
        assert_almost_equal(54.7170296089477, intersection_point_2.lat().v())
        assert_almost_equal(-14.56385574430775, intersection_point_2.lon().v())


def test_intersection_point_distance_non_wgs84():
    # Example from Charles Karney email on 31/03/2025
    # Issue #1 Python calculate_intersection_point code crashes
    ell = Ellipsoid(Metres(6.4e6), 1.0/50.0)
    g1 = Geodesic(LatLong(Degrees(-30), Degrees(0.0)),
                  LatLong(Degrees(29.5), Degrees(179.5)),
                  ell)
    print(g1.length())
    g2 = Geodesic(LatLong(Degrees(1), Degrees(90)),
                  LatLong(Degrees(-2), Degrees(-95)),
                  ell)
    print(g2.length())
    p = calculate_intersection_point(g1, g2, Metres(1e-6))
    print(p)

if __name__ == '__main__':
    pytest.main()
