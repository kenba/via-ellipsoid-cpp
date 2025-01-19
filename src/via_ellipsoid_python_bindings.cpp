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
//
/// @file via_ellipsoid_python_bindings.cpp
/// @brief Contains the via::ellipsoid python interface
//////////////////////////////////////////////////////////////////////////////
// ensure numpy.h included before ellipsoid.hpp
// clang-format off
#include <pybind11/numpy.h>
#include "via/ellipsoid.hpp"
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
      .def("a3", &EllipsoidDouble::a3)
      .def("c3x", &EllipsoidDouble::c3x)
      .def("calculate_epsilon", &EllipsoidDouble::calculate_epsilon)
      .def("calculate_parametric_latitude",
           &EllipsoidDouble::calculate_parametric_latitude)
      .def("calculate_geodetic_latitude",
           &EllipsoidDouble::calculate_geodetic_latitude);

  // Python bindings for geodesic functions
  m.def("find_azimuth_and_aux_length",
        &via::ellipsoid::find_azimuth_and_aux_length<double>,
        "Calculate the initial azimuth and great circle length between a pair "
        "of points on the auxiliary sphere.");
  m.def("aux_sphere_azimuth_length",
        &via::ellipsoid::aux_sphere_azimuth_length<double>,
        "Calculate the initial azimuth and great circle length between a pair "
        "of points on the auxiliary sphere.");
  m.def("calculate_azimuth_aux_length",
        &via::ellipsoid::calculate_azimuth_aux_length<double>,
        "Calculate the `geodesic` azimuth and great circle length on the "
        "auxiliary sphere between a pair of positions.");
  m.def("convert_radians_to_metres",
        &via::ellipsoid::convert_radians_to_metres<double>,
        "Convert a great circle distance on the auxiliary sphere in radians to "
        "metres on the ellipsoid.");

  // Python numpy binding for the Geodesic class
  using GeodesicDouble = via::ellipsoid::Geodesic<double>;
  py::class_<GeodesicDouble>(m, "Geodesic")
      .def(py::init<via::LatLong<double>, via::Angle<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::units::si::Metres<double>>())
      .def(py::init<via::LatLong<double>, via::Angle<double>,
                    via::units::si::Metres<double>, EllipsoidDouble>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>>())
      .def(py::init<via::LatLong<double>, via::LatLong<double>,
                    EllipsoidDouble>())
      .def("is_valid", &GeodesicDouble::is_valid)
      .def("beta", &GeodesicDouble::beta)
      .def("lon", &GeodesicDouble::lon)
      .def("azi", &GeodesicDouble::azi)
      .def("set_aux_length", &GeodesicDouble::set_aux_length)
      .def("aux_length", &GeodesicDouble::aux_length)
      .def("ellipsoid", &GeodesicDouble::ellipsoid)
      .def("a", &GeodesicDouble::a)
      .def("pole", &GeodesicDouble::pole)
      .def("direction", &GeodesicDouble::direction)
      .def("clairaut", &GeodesicDouble::clairaut)
      .def("epsilon", &GeodesicDouble::epsilon)
      .def("metres_to_radians", &GeodesicDouble::metres_to_radians)
      .def("radians_to_metres", &GeodesicDouble::radians_to_metres)
      .def("length", &GeodesicDouble::length)
      .def("aux_beta", &GeodesicDouble::aux_beta)
      .def("aux_latitude", &GeodesicDouble::aux_latitude)
      .def("latitude", &GeodesicDouble::latitude)
      .def("aux_azimuth", &GeodesicDouble::aux_azimuth)
      .def("azimuth", &GeodesicDouble::azimuth)
      .def("delta_longitude", &GeodesicDouble::delta_longitude)
      .def("aux_longitude", &GeodesicDouble::aux_longitude)
      .def("longitude", &GeodesicDouble::longitude)
      .def("aux_lat_long", &GeodesicDouble::aux_lat_long)
      .def("lat_long", &GeodesicDouble::lat_long)
      .def("aux_point", &GeodesicDouble::aux_point)
      .def("b", &GeodesicDouble::b)
      .def("mid_point", &GeodesicDouble::mid_point)
      .def("aux_point_and_pole", &GeodesicDouble::aux_point_and_pole)
      .def("calculate_aux_atd_and_xtd",
           &GeodesicDouble::calculate_aux_atd_and_xtd<10>)
      .def("calculate_atd_and_xtd", &GeodesicDouble::calculate_atd_and_xtd<10>);

  // Python bindings for Geodesic intersection functions
  m.def("calculate_geodesic_intersection_distances",
        &via::ellipsoid::calculate_geodesic_intersection_distances<double>,
        "Calculate the auxiliary Great Circle arc lengths to an intersection "
        "point of two geodesics.");
  m.def("calculate_aux_intersection_distances",
        &via::ellipsoid::calculate_aux_intersection_distances<double>,
        "Calculate the distances along a pair of Geodesics (in Radians) to "
        "their closest intersection or reference points.");
  m.def("calculate_intersection_distances",
        &via::ellipsoid::calculate_intersection_distances<double>,
        "Calculate the distances along a pair of Geodesics (in Radians) to "
        "their closest intersection or reference points.");
  m.def("calculate_intersection_point",
        &via::ellipsoid::calculate_intersection_point<double>,
        "Calculate the position (Latitude and Longitude) where a pair of "
        "`Geodesic`s intersect, or None if the `Geodesic`s do not intersect.");
}
