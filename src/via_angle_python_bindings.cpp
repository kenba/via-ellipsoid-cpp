//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2024 Ken Barker
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
/// @file via_angle_python_bindings.cpp
/// @brief Contains the via::angle python interface
//////////////////////////////////////////////////////////////////////////////
// ensure numpy.h included before angle.hpp
// clang-format off
#include <pybind11/numpy.h>
#include "via/angle.hpp"
// clang-format on
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(via_angle, m) {
  try {
    py::module::import("numpy");
  } catch (...) {
    return;
  }

  // Python bindings for constants
  m.attr("SQ_EPSILON") = via::trig::SQ_EPSILON<double>;
  m.attr("TAU") = via::trig::TAU<double>;
  m.attr("PI_3") = via::trig::PI_3<double>;
  m.attr("SQRT1_2") = via::trig::SQRT1_2<double>;
  m.attr("COS_30_DEGREES") = via::trig::COS_30_DEGREES<double>;

  // Python bindings for degrees /radians conversion functions
  m.def("deg2rad", &via::trig::deg2rad<double>,
        "Convert a value in degrees to radians.");
  m.def("rad2deg", &via::trig::rad2deg<double>,
        "Convert a value in radians to degrees.");

  // Python numpy binding for the UnitNegRange class
  PYBIND11_NUMPY_DTYPE(via::trig::UnitNegRange<double>, v_);

  // Python bindings for the UnitNegRange class
  py::class_<via::trig::UnitNegRange<double>>(m, "UnitNegRange")
      .def(py::init<double>())

      .def("v", &via::trig::UnitNegRange<double>::v)
      .def("clamp", &via::trig::UnitNegRange<double>::clamp)
      .def("__abs__", &via::trig::UnitNegRange<double>::abs)
      .def("__repr__", &via::trig::UnitNegRange<double>::python_repr)

      .def(-py::self)
      .def(py::self == py::self);

  // Python bindings for the trig functions
  m.def("swap_sin_cos", &via::trig::swap_sin_cos<double>,
        "Convert a cosine to a sine, or vice versa.");
  m.def("cosine_from_sine", &via::trig::cosine_from_sine<double>,
        "Convert a cosine to a sine, or vice versa.");
  m.def("sine", &via::trig::sine<double>,
        "Calculate the sine of an angle in radians.");
  m.def("cosine", &via::trig::cosine<double>,
        "Calculate the cosine of an angle in radians using the sine of the "
        "angle.");
  m.def("sincos", &via::trig::sincos<double>,
        "Calculate the sine and cosine of an angle from a value in radians.");
  m.def("sincos_diff", &via::trig::sincos_diff<double>,
        "Calculate the sine and cosine of an angle from a value in radians.");
  m.def("arctan2", &via::trig::arctan2<double>,
        "Calculate an angle in radians from its sine and cosine.");
  m.def("sincosd", &via::trig::sincosd<double>,
        "Calculate the sine and cosine of an angle from a value in degrees.");
  m.def("sincosd_diff", &via::trig::sincosd_diff<double>,
        "Calculate the sine and cosine of an angle from the difference of a "
        "pair of values in degrees.");
  m.def("arctan2d", &via::trig::arctan2d<double>,
        "Calculate an angle in degrees from its sine and cosine.");

  m.def("calculate_adjacent_length",
        &via::trig::calculate_adjacent_length<double>,
        "Calculates the length of the other side in a right angled triangle.");

  m.def("spherical_adjacent_length",
        &via::trig::spherical_adjacent_length<double>,
        "Calculates the length of the other side in a right angled spherical "
        "triangle.");

  m.def("spherical_hypotenuse_length",
        &via::trig::spherical_hypotenuse_length<double>,
        "Calculates the length of the hypotenuse in a right angled spherical "
        "triangle.");

  m.def(
      "spherical_cosine_rule", &via::trig::spherical_cosine_rule<double>,
      "Calculates the length of the adjacent side of a right angled spherical "
      "triangle.");

  // Python numpy binding for the Degrees class
  PYBIND11_NUMPY_DTYPE(via::Degrees<double>, v_);

  // Python bindings for the Degrees class
  py::class_<via::Degrees<double>>(m, "Degrees")
      .def(py::init<double>())

      .def("v", &via::Degrees<double>::v)
      .def("opposite", &via::Degrees<double>::opposite)
      .def("__abs__", &via::Degrees<double>::abs)
      .def("__repr__", &via::Degrees<double>::python_repr)

      .def(-py::self)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self == py::self);

  // Python numpy binding for the Radians class
  PYBIND11_NUMPY_DTYPE(via::Radians<double>, v_);

  // Python bindings for the Radians class
  py::class_<via::Radians<double>>(m, "Radians")
      .def(py::init<double>())

      .def("v", &via::Radians<double>::v)
      .def("clamp", &via::Radians<double>::opposite)
      .def("opposite", &via::Radians<double>::opposite)
      .def("__abs__", &via::Radians<double>::abs)
      .def("__repr__", &via::Radians<double>::python_repr)

      .def(-py::self)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self == py::self);

  // Python numpy binding for the Angle class
  PYBIND11_NUMPY_DTYPE(via::Angle<double>, sin_, cos_);

  // Python bindings for the Angle class
  py::class_<via::Angle<double>>(m, "Angle")
      .def(py::init<>())
      .def(py::init<via::trig::UnitNegRange<double>,
                    via::trig::UnitNegRange<double>>())
      .def(py::init<via::Degrees<double>>())
      .def(py::init<via::Degrees<double>, via::Degrees<double>>())
      .def(py::init<via::Radians<double>>())
      .def(py::init<via::Radians<double>, via::Radians<double>>())

      .def("from_y_x", via::Angle<double>::from_y_x)

      .def("is_valid", &via::Angle<double>::is_valid)

      .def("sin", &via::Angle<double>::sin)
      .def("cos", &via::Angle<double>::cos)
      .def("tan", &via::Angle<double>::tan)
      .def("csc", &via::Angle<double>::csc)
      .def("sec", &via::Angle<double>::sec)
      .def("cot", &via::Angle<double>::cot)

      .def("to_degrees", &via::Angle<double>::to_degrees)
      .def("to_radians", &via::Angle<double>::to_radians)
      .def("opposite", &via::Angle<double>::opposite)
      .def("quarter_turn_cw", &via::Angle<double>::quarter_turn_cw)
      .def("quarter_turn_ccw", &via::Angle<double>::quarter_turn_ccw)
      .def("negate_cos", &via::Angle<double>::negate_cos)
      .def("x2", &via::Angle<double>::x2)
      .def("half", &via::Angle<double>::half)

      .def("__abs__", &via::Angle<double>::abs)
      .def("__repr__", &via::Angle<double>::python_repr)

      .def(-py::self)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self == py::self);

  m.def("sine_sum", &via::sine_sum<double>,
        "Calculate the sine of the sum of two angles.");
  m.def("sine_diff", &via::sine_diff<double>,
        "Calculate the sine of the difference of two angles.");
  m.def("cosine_sum", &via::cosine_sum<double>,
        "Calculate the cosine of the sum of two angles.");
  m.def("cosine_diff", &via::cosine_diff<double>,
        "Calculate the cosine of the difference of two angles.");
}
