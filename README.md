# via-ellipsoid-cpp

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/license/mit/)
[![C/C++ CI](https://github.com/kenba/via-ellipsoid-cpp/workflows/C/C++%20CI/badge.svg)](https://github.com/kenba/via-ellipsoid-cpp/actions)
[![codecov](https://codecov.io/gh/kenba/via-ellipsoid-cpp/graph/badge.svg?token=CIPRUEW7RT)](https://codecov.io/gh/kenba/via-ellipsoid-cpp)

A library for performing geometric calculations on the
[WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System) ellipsoid,
see *Figure 1*.

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/3e/WGS84_mean_Earth_radius.svg/800px-WGS84_mean_Earth_radius.svg.png" width="400">

*Figure 1 The WGS84 Ellipsoid (not to scale)*

WGS84 has become the de facto standard for satellite navigation since its adoption
by the Navstar [Global Positioning System](https://en.wikipedia.org/wiki/Global_Positioning_System)
(GPS) and US president Ronald Reagan's 1983 decision to make GPS available for civilian use
after airliner [KAL 007](https://en.wikipedia.org/wiki/Korean_Air_Lines_Flight_007)
was shot down by Soviet interceptor aircraft when it strayed into
prohibited airspace due to navigational errors.

This library uses the WGS84 primary parameters defined in Tab. 3-1 of the
[ICAO WGS 84 Implementation Manual](https://www.icao.int/safety/pbn/Documentation/EUROCONTROL/Eurocontrol%20WGS%2084%20Implementation%20Manual.pdf).

## Geodesic navigation

The shortest path between two points on the surface of an ellipsoid is a
[geodesic](https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid) -
the equivalent of straight line segments in planar geometry or
[great circles](https://en.wikipedia.org/wiki/Great_circle) on the surface of a
sphere, see *Figure 2*.

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Geodesic_problem_on_an_ellipsoid.svg/1024px-Geodesic_problem_on_an_ellipsoid.svg.png" width="400">

*Figure 2 A geodesic between points A and B*

This library uses the correspondence between geodesics on an ellipsoid
and great-circles on the auxiliary sphere together with 3D vectors to calculate:

- the initial azimuth and length of a geodesic between two positions;
- the along track distance and across track distance of a position relative to a geodesic;
- and the intersection of a pair of geodesics.
