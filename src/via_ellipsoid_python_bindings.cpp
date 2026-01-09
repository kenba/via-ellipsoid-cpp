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
//
/// @file via_ellipsoid_python_bindings.cpp
/// @brief Contains the via::ellipsoid python interface
//////////////////////////////////////////////////////////////////////////////
// ensure numpy.h included before ellipsoid.hpp
// clang-format off
#include <pybind11/numpy.h>
#include "via/ellipsoid.hpp"
#include "via/ellipsoid/vincenty_functions.hpp"
// clang-format on
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(via_ellipsoid, m) {
  try {
    py::module::import("numpy");
  } catch (...) {
    return;
  }

  // Python bindings for constants
  m.attr("WGS84_A") = via::ellipsoid::wgs84::A<double>;
  m.attr("WGS84_F") = via::ellipsoid::wgs84::F<double>;

  // Python numpy binding for the Ellipsoid class
  using EllipsoidDouble = via::ellipsoid::Ellipsoid<double>;
  PYBIND11_NUMPY_DTYPE(EllipsoidDouble, a_, f_, b_, one_minus_f_,
                       recip_one_minus_f_, e_2_, ep_2_, n_, a3_, c3x_);

  // Python bindings for the Ellipsoid class
  py::class_<EllipsoidDouble>(m, "Ellipsoid")
      .def_static("wgs84", &EllipsoidDouble::wgs84)
      .def(py::init<via::units::si::Metres<double>, double>())
      .def("a", &EllipsoidDouble::a)
      .def("f", &EllipsoidDouble::f)
      .def("b", &EllipsoidDouble::b)
      .def("one_minus_f", &EllipsoidDouble::one_minus_f)
      .def("recip_one_minus_f", &EllipsoidDouble::recip_one_minus_f)
      .def("e_2", &EllipsoidDouble::e_2)
      .def("ep_2", &EllipsoidDouble::ep_2)
      .def("n", &EllipsoidDouble::n)
      .def("calculate_epsilon", &EllipsoidDouble::calculate_epsilon)
      .def("calculate_a3f", &EllipsoidDouble::calculate_a3f)
      .def("calculate_a3c", &EllipsoidDouble::calculate_a3c)
      .def("calculate_c3y", &EllipsoidDouble::calculate_c3y)
      .def("calculate_parametric_latitude",
           &EllipsoidDouble::calculate_parametric_latitude)
      .def("calculate_geodetic_latitude",
           &EllipsoidDouble::calculate_geodetic_latitude);

  // Python bindings for geodesic functions
  m.def("find_azimuths_and_arc_length",
        &via::ellipsoid::find_azimuths_and_arc_length<double>,
        "Find the aziumth and arc length on the auxiliary sphere.");
  m.def("aux_sphere_azimuths_length",
        &via::ellipsoid::aux_sphere_azimuths_length<double>,
        "Calculate the initial azimuth and great circle length between a pair "
        "of points on the auxiliary sphere.");
  m.def("calculate_azimuths_arc_length",
        &via::ellipsoid::calculate_azimuths_arc_length<double>,
        "Calculate the `geodesic` azimuth and great circle length on the "
        "auxiliary sphere between a pair of positions.");
  m.def("convert_radians_to_metres",
        &via::ellipsoid::convert_radians_to_metres<double>,
        "Convert a great circle distance on the auxiliary sphere in radians to "
        "metres on the ellipsoid.");

  // Python bindings for geodesic vincenty functions
  m.def("inverse_azimuths_and_distance",
        &via::ellipsoid::vincenty::inverse_azimuths_and_distance<double>,
        "Calculate the `geodesic` azimuths and distance between a pair of "
        "positions using Vincenty's algorithm.");

  // Python numpy binding for the GeodesicSegment class
  using GeodesicSegmentDouble = via::ellipsoid::GeodesicSegment<double>;
  py::class_<GeodesicSegmentDouble>(m, "GeodesicSegment")
      .def(py::init<via::LatLong<double>, via::Angle<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::Radians<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::Radians<double>, via::units::si::Metres<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::Radians<double>, via::units::si::Metres<double>,
                    const EllipsoidDouble &>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::units::si::Metres<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::units::si::Metres<double>,
                    via::units::si::Metres<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::units::si::Metres<double>,
                    via::units::si::Metres<double>, const EllipsoidDouble &>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>,
                    via::units::si::Metres<double>>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>,
                    via::units::si::Metres<double>, via::Radians<double>>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>,
                    via::units::si::Metres<double>, via::Radians<double>,
                    const EllipsoidDouble &>())
      .def("is_valid", &GeodesicSegmentDouble::is_valid)
      .def("beta", &GeodesicSegmentDouble::beta)
      .def("lon", &GeodesicSegmentDouble::lon)
      .def("azi", &GeodesicSegmentDouble::azi)
      .def("set_arc_length", &GeodesicSegmentDouble::set_arc_length)
      .def("arc_length", &GeodesicSegmentDouble::arc_length)
      .def("set_half_width", &GeodesicSegmentDouble::set_half_width)
      .def("half_width", &GeodesicSegmentDouble::half_width)
      .def("ellipsoid", &GeodesicSegmentDouble::ellipsoid)
      .def("epsilon", &GeodesicSegmentDouble::epsilon)
      .def("a", &GeodesicSegmentDouble::a)
      .def("metres_to_radians", &GeodesicSegmentDouble::metres_to_radians)
      .def("length", &GeodesicSegmentDouble::length)
      .def("arc_beta", &GeodesicSegmentDouble::arc_beta)
      .def("arc_latitude", &GeodesicSegmentDouble::arc_latitude)
      .def("latitude", &GeodesicSegmentDouble::latitude)
      .def("arc_azimuth", &GeodesicSegmentDouble::arc_azimuth)
      .def("azimuth", &GeodesicSegmentDouble::azimuth)
      .def("delta_longitude", &GeodesicSegmentDouble::delta_longitude)
      .def("arc_longitude", &GeodesicSegmentDouble::arc_longitude)
      .def("longitude", &GeodesicSegmentDouble::longitude)
      .def("arc_beta_long", &GeodesicSegmentDouble::arc_beta_long)
      .def("arc_lat_long", &GeodesicSegmentDouble::arc_lat_long)
      .def("lat_long", &GeodesicSegmentDouble::lat_long)
      .def("arc_angles", &GeodesicSegmentDouble::arc_angles)
      .def("arc_point", &GeodesicSegmentDouble::arc_point)
      .def("mid_point", &GeodesicSegmentDouble::mid_point)
      .def("arc_pole", &GeodesicSegmentDouble::arc_pole)
      .def("arc_point_and_pole", &GeodesicSegmentDouble::arc_point_and_pole)
      .def("calculate_arc_atd_and_xtd",
           &GeodesicSegmentDouble::calculate_arc_atd_and_xtd<10>)
      .def("calculate_sphere_shortest_distance",
           &GeodesicSegmentDouble::calculate_sphere_shortest_distance<10>)
      .def("calculate_atd_and_xtd",
           &GeodesicSegmentDouble::calculate_atd_and_xtd<10>)
      .def("shortest_distance", &GeodesicSegmentDouble::shortest_distance<10>)

      .def(py::self == py::self);

  // Python bindings for GeodesicSegment intersection functions
  m.def("geodesics_are_coincident",
        &via::ellipsoid::intersection::geodesics_are_coincident<double>,
        "Determine whether two `GeodesicSegment`s are coincident.");
  m.def("calculate_arc_reference_distances_and_angle",
        &via::ellipsoid::intersection::
            calculate_arc_reference_distances_and_angle<double>,
        "Find the closest intersection distances of two `GeodesicSegment`s and "
        "the relative angle at the reference point.");

  m.def("calculate_intersection_distances",
        &via::ellipsoid::calculate_intersection_distances<double>,
        "Calculate the distances along a pair of GeodesicSegments (in Radians) "
        "to their closest intersection or reference points.");
  m.def("calculate_intersection_point",
        &via::ellipsoid::calculate_intersection_point<double>,
        "Calculate the position (Latitude and Longitude) where a pair of "
        "`GeodesicSegment`s intersect, or None if the `GeodesicSegment`s do "
        "not intersect.");
}
