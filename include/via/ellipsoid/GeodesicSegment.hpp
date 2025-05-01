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
/// @file GeodesicSegment.hpp
/// @brief Contains the via::ellipsoid GeodesicSegment class.
//////////////////////////////////////////////////////////////////////////////
#include "geodesic_functions.hpp"
#include <via/angle/trig.hpp>
#include <via/sphere/great_circle.hpp>

namespace via {
namespace ellipsoid {
/// A geodesic segment on the surface of an ellipsoid.
///
/// A geodesic segment on an ellipsoid is the shortest path between two points.
/// It is represented by a great circle arc on the auxiliary sphere.
/// @invariant -90° <= beta_ <= 90°
/// @invariant 0 <= arc_length_ <= π
template <typename T, typename V = vector::Vector3<T>>
  requires std::floating_point<T>
class GeodesicSegment final {
#ifdef PYBIND11_NUMPY_DTYPE
public:
#endif
  /// The parametric start latitude on the auxiliary sphere.
  Angle<T> beta_ = Angle<T>(T(-1), T());
  Angle<T> lon_ = Angle<T>(T(), T(1)); ///< The start longitude.
  Angle<T> azi_ = Angle<T>(T(), T(1)); ///< The start azimuth.

  /// Azimuth at the Equator.
  Angle<T> azi0_ = Angle<T>(T(), T(1));
  /// Great Circle distance from Northward Equator crossing.
  Angle<T> sigma1_ = Angle<T>(T(), T(1));
  /// The Great Circle length on the auxiliary sphere in radians.
  Radians<T> arc_length_ = Radians<T>(T());
  /// integration constant epsilon, derived from Clairaut's constant.
  T eps_ = T();
  T a1_ = T(1); ///< constant used to convert geodesic/great circle distances.
  T a3c_ = T(); ///< constant used to convert geodesic/great circle longitudes.
  /// start point geodesic/great circle distance difference.
  Radians<T> b11_ = Radians<T>(T());
  /// A const reference to the underlying Ellipsoid, default WGS84.
  const Ellipsoid<T> &ellipsoid_{Ellipsoid<T>::wgs84()};

#ifndef PYBIND11_NUMPY_DTYPE
public:
#endif

  /// Construct a GeodesicSegment from a start point, parametric lat lon,
  /// azimuth and great circle arc length on the auxiliary sphere.
  /// @pre -90° <= beta <= 90°
  /// @pre 0 <= arc_length <= π
  /// @param beta the parametric start latitude.
  /// @param lon the start longitude.
  /// @param azimuth the azimuth at the start point.
  /// @param arc_length the great circle arc length in Radians.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const Angle<T> beta, const Angle<T> lon, const Angle<T> azimuth,
      const Radians<T> arc_length,
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : beta_{beta}, lon_{lon}, azi_{azimuth},
        // Calculate the azimuth at the first Equator crossing
        azi0_(trig::UnitNegRange(azi_.sin().v() * beta_.cos().v()),
              trig::swap_sin_cos(
                  trig::UnitNegRange(azi_.sin().v() * beta_.cos().v()))),
        // Calculate the distance to the first Equator crossing
        sigma1_{Angle<T>::from_y_x(beta_.sin().v(),
                                   beta.cos().v() * azi_.cos().v())},
        arc_length_{arc_length},
        // Calculate eps for calculating coefficients
        eps_(ellipsoid.calculate_epsilon(azi0_.sin())),
        a1_((evaluate_A1<T>(eps_) + T(1))),
        a3c_(ellipsoid.f() * azi0_.sin().v() *
             evaluate_poynomial(eps_, ellipsoid.a3())),
        ellipsoid_(ellipsoid) {
    Expects((T() <= beta.cos().v()) && (T() <= arc_length.v()) &&
            (arc_length.v() <= trig::PI<T>));

    const auto C1{evaluate_coeffs_C1<T>(eps_)};
    b11_ = sin_cos_series(sigma1_, C1);
  }

  /// Construct a GeodesicSegment from a start point, azimuth and great circle
  /// arc length on the auxiliary sphere.
  /// @pre 0 <= arc_length <= π
  /// @param a the geodetic start point.
  /// @param azimuth the azimuth at the start point.
  /// @param arc_length the great circle arc length in Radians.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const LatLong<T> &a, const Angle<T> &azimuth,
      const Radians<T> arc_length = Radians(T()),
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : GeodesicSegment(
            ellipsoid.calculate_parametric_latitude(Angle<T>(a.lat())),
            Angle<T>(a.lon()), azimuth, arc_length, ellipsoid) {}

  /// Construct a GeodesicSegment from a start point, azimuth and length in
  /// Metres.
  /// @param a the geodetic start point.
  /// @param azimuth the azimuth at the start point.
  /// @param length the length in Metres.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const LatLong<T> &a, const Angle<T> &azimuth,
      const units::si::Metres<T> length,
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : GeodesicSegment(
            ellipsoid.calculate_parametric_latitude(Angle<T>(a.lat())),
            Angle<T>(a.lon()), azimuth, Radians(T()), ellipsoid) {
    arc_length_ = metres_to_radians(length);
  }

  /// Construct a GeodesicSegment from a start point and azimuth, arc length
  /// tuple.
  /// @param a the geodetic point.
  /// @param azimuth_length the azimuth, arc length tuple.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const LatLong<T> &a,
      const std::tuple<Angle<T>, Radians<T>, Angle<T>, unsigned> azimuth_length,
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : GeodesicSegment(a, std::get<0>(azimuth_length),
                        std::get<1>(azimuth_length), ellipsoid) {}

  /// Construct a GeodesicSegment between a pair of points.
  ///
  /// @pre a.is_valid() && b.s_valid()
  /// @pre tolerance >= epsilon
  ///
  /// @param a, b the start and finish points in geodetic coordinates.
  /// @param tolerance the tolerance to perform the calculation to in Radians,
  /// default great_circle::MIN_VALUE.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const LatLong<T> &a, const LatLong<T> &b,
      const Radians<T> tolerance = Radians<T>(great_circle::MIN_VALUE<T>),
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : GeodesicSegment(
            a, calculate_azimuths_arc_length(a, b, tolerance, ellipsoid),
            ellipsoid) {}

  /// Test whether a `GeodesicSegment` is valid:
  /// @invariant 0° <= latitude <= 90°
  /// @invariant 0 <= arc_length <= π
  [[nodiscard("Pure Function")]]
  constexpr auto is_valid() const noexcept -> bool {
    return beta_.cos().v() >= T() && T() <= arc_length_.v() &&
           arc_length_.v() <= trig::PI<T>;
  }

  /// Accessor for the start parametric latitude on the auxiliary sphere.
  /// @return beta_.
  [[nodiscard("Pure Function")]]
  constexpr auto beta() const noexcept -> Angle<T> {
    return beta_;
  }

  /// Accessor for the start longitude.
  /// @return lon_.
  [[nodiscard("Pure Function")]]
  constexpr auto lon() const noexcept -> Angle<T> {
    return lon_;
  }

  /// Accessor for the start azimuth.
  /// @return azi_.
  [[nodiscard("Pure Function")]]
  constexpr auto azi() const noexcept -> Angle<T> {
    return azi_;
  }

  /// Accessor for the arc length on the auxiliary sphere in radians.
  /// @return The arc length on the auxiliary sphere in radians.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_length() const noexcept -> Radians<T> {
    return arc_length_;
  }

  /// Accessor for the reference to the underlying `Ellipsoid`.
  /// @return a const reference to the Ellipsoid.
  [[nodiscard("Pure Function")]]
  constexpr auto ellipsoid() const noexcept -> const Ellipsoid<T> & {
    return ellipsoid_;
  }

  /// Accessor for the integration constant: epsilon.
  /// @return eps_.
  [[nodiscard("Pure Function")]]
  constexpr auto epsilon() const noexcept -> T {
    return eps_;
  }

  /// Accessor for the start point on the unit sphere.
  /// @return start point
  [[nodiscard("Pure Function")]]
  constexpr auto a() const noexcept -> V {
    return vector::to_point<T>(beta_, lon_);
  }

  /// Convert a distance in metres on the ellipsoid to radians on the
  /// auxiliary sphere.
  /// @param distance_m the distance along the GeodesicSegment in metres.
  /// @return the distance along the great circle arc in radians.
  [[nodiscard("Pure Function")]]
  constexpr auto metres_to_radians(const units::si::Metres<T> distance_m) const
      -> Radians<T> {
    if (std::abs(distance_m.v()) < great_circle::MIN_VALUE<T>)
      return Radians(T());

    const auto tau12{Radians<T>(distance_m.v() / (ellipsoid_.b().v() * a1_))};
    const auto tau_sum{Angle<T>(b11_ + tau12)};
    const auto c1p{evaluate_coeffs_C1p<T>(eps_)};
    const auto b12{sin_cos_series(sigma1_ + Angle<T>(tau_sum), c1p)};

    return tau12 + b12 + b11_;
  }

  /// Convert a distance in radians on the auxiliary sphere to metres
  /// on the ellipsoid.
  /// @param gc_distance the distance along the great circle arc in radians.
  /// @return the distance along the GeodesicSegment in metres.
  [[nodiscard("Pure Function")]]
  constexpr auto radians_to_metres(const Radians<T> gc_distance) const
      -> units::si::Metres<T> {
    if (gc_distance.abs().v() < great_circle::MIN_VALUE<T>)
      return units::si::Metres(T());

    // Calculate Great circle distance from Northward Equator crossing.
    const Angle<T> sigma_sum{sigma1_ + Angle<T>(gc_distance)};
    const auto c1{evaluate_coeffs_C1<T>(eps_)};
    const auto b12{sin_cos_series(sigma_sum, c1)};

    return units::si::Metres<T>(ellipsoid().b().v() * a1_ *
                                (gc_distance + b12 - b11_).v());
  }

  /// Accessor for the length of the GeodesicSegment in metres.
  /// @return The length of the GeodesicSegment in metres.
  [[nodiscard("Pure Function")]]
  constexpr auto length() const -> units::si::Metres<T> {
    return radians_to_metres(arc_length_);
  }

  /// Calculate the parametric latitude at the great circle length.
  /// @param length the length on the auxiliary sphere as an Angle.
  /// @return the parametric latitude.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_beta(const Angle<T> length) const -> Angle<T> {
    const trig::UnitNegRange sin_beta{length.cos().v() * beta_.sin().v() +
                                      length.sin().v() * beta_.cos().v() *
                                          azi_.cos().v()};
    const Angle<T> beta(sin_beta, trig::swap_sin_cos(sin_beta));

    Ensures(beta.is_valid());

    return beta;
  }

  /// Calculate the geodetic latitude at the great circle length.
  /// @param gc_length the great circle length on the auxiliary sphere, in
  /// radians.
  /// @return the geodetic latitude of the position at `gc_length`.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_latitude(const Radians<T> gc_length) const -> Angle<T> {
    return ellipsoid_.calculate_geodetic_latitude(
        arc_beta(Angle<T>(gc_length)));
  }

  /// Calculate the geodetic latitude at the length along the geodesic.
  /// @param length_m the length along the geodesic, in metres.
  /// @return the geodetic latitude.
  [[nodiscard("Pure Function")]]
  constexpr auto latitude(const units::si::Metres<T> length_m) const
      -> Angle<T> {
    if (std::abs(length_m.v()) < great_circle::MIN_VALUE<T>)
      return ellipsoid_.calculate_geodetic_latitude(beta_);

    return arc_latitude(metres_to_radians(length_m));
  }

  /// Calculate the azimuth along the arc at great circle distance sigma.
  /// @param sigma the great circle distance on the auxiliary sphere.
  /// @return the azimuth of the geodesic/great circle at `gc_length`.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_azimuth(const Angle<T> sigma) const -> Angle<T> {
    constexpr T MAX_LAT{T(1) - std::numeric_limits<T>::epsilon()};

    // Use Karney's method to calculate latitude and azimuth.
    const Angle<T> sigma_sum{sigma1_ + sigma};
    const T sin_beta{azi0_.cos().v() * sigma_sum.sin().v()};

    // Handle North pole, only valid azimuth is South
    if (sin_beta > MAX_LAT)
      return Angle<T>(trig::UnitNegRange(T()), trig::UnitNegRange(T(-1)));
    else
      return Angle<T>::from_y_x(azi0_.sin().v(),
                                azi0_.cos().v() * sigma_sum.cos().v());
  }

  /// Calculate the azimuth at the length along the geodesic.
  /// @param length_m the length along the geodesic segment, in metres.
  /// @return the azimuth of the geodesic segment at length.
  [[nodiscard("Pure Function")]]
  constexpr auto
  azimuth(const units::si::Metres<T> length_m = units::si::Metres(T())) const
      -> Angle<T> {
    const Angle<T> sigma{metres_to_radians(length_m)};
    return arc_azimuth(sigma);
  }

  /// Calculate the geodesic longitude difference at a great circle length
  /// along the auxiliary sphere.
  /// @param gc_length the great circle length on the auxiliary sphere,
  /// in radians.
  /// @return the longitude difference from the start point.
  [[nodiscard("Pure Function")]]
  constexpr auto delta_longitude(const Radians<T> gc_length) const -> Angle<T> {
    if (gc_length.abs().v() <= great_circle::MIN_VALUE<T>)
      return Angle<T>();

    // The great circle distance from Northward Equator crossing.
    const Angle<T> sigma_sum(sigma1_ + Angle<T>(gc_length));

    // The longitude difference on the auxiliary sphere.
    const auto omega1{Angle<T>::from_y_x(azi0_.sin().v() * beta_.sin().v(),
                                         beta_.cos().v() * azi_.cos().v())};
    const auto omega2{Angle<T>::from_y_x(azi0_.sin().v() * sigma_sum.sin().v(),
                                         sigma_sum.cos().v())};
    const Angle<T> omega12{omega2 - omega1};

    const auto c3{evaluate_coeffs_C3y<T>(ellipsoid_.c3x(), eps_)};
    const auto b31(sin_cos_series(sigma1_, c3));
    const auto b32(sin_cos_series(sigma_sum, c3));

    return omega12 -
           Angle<T>(Radians<T>{a3c_ * (gc_length.v() + (b32.v() - b31.v()))});
  }

  /// Calculate the geodesic longitude at the great circle length along
  /// the auxiliary sphere.
  /// @param gc_length the great circle length on the auxiliary sphere,
  /// in radians.
  /// @return the geodesic longitude as an Angle.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_longitude(const Radians<T> gc_length) const -> Angle<T> {
    return lon_ + delta_longitude(gc_length);
  }

  /// Calculate the geodesic longitude at the length along the geodesic.
  /// @param length_m the length along the geodesic, in metres.
  /// @return the geodesic longitude as an Angle.
  [[nodiscard("Pure Function")]]
  constexpr auto longitude(const units::si::Metres<T> length_m) const
      -> Angle<T> {
    return arc_longitude(metres_to_radians(length_m));
  }

  /// Calculate the geodesic `LatLong` at the great circle length along
  /// the auxiliary sphere.
  /// @param gc_length the great circle length on the auxiliary sphere,
  /// in Radians.
  /// @return the `LatLong` of the geodesic position at `gc_length`.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_lat_long(const Radians<T> gc_length) const -> LatLong<T> {
    return LatLong<T>(arc_latitude(gc_length).to_degrees(),
                      arc_longitude(gc_length).to_degrees());
  }

  /// Calculate the geodesic `LatLong` at the length along the
  /// `GeodesicSegment`.
  /// @param length_m the length along the geodesic, in metres.
  /// @return the `LatLong` of the geodesic position at `length_m`.
  [[nodiscard("Pure Function")]]
  constexpr auto lat_long(const units::si::Metres<T> length_m) const
      -> LatLong<T> {
    return arc_lat_long(metres_to_radians(length_m));
  }

  /// Calculate the point on the auxiliary sphere at the
  /// great circle arc length.
  /// @param gc_length the great circle arc length, in radians.
  /// @return the point on the auxiliary sphere at `gc_length`.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_point(const Radians<T> gc_length) const -> V {
    if (gc_length.abs().v() < great_circle::MIN_VALUE<T>)
      return a();

    const Angle<T> beta{arc_beta(Angle<T>(gc_length))};
    const Angle<T> lon{arc_longitude(gc_length)};
    return vector::to_point(beta, lon);
  }

  /// Calculate the end point on the geodesic.
  /// @return end point
  [[nodiscard("Pure Function")]]
  constexpr auto b() const -> V {
    return arc_point(arc_length_);
  }

  /// Calculate the point on the auxiliary sphere at the mid point of the
  /// `GeodesicSegment`.
  /// @return end point
  [[nodiscard("Pure Function")]]
  constexpr auto mid_point() const -> V {
    return arc_point(
        metres_to_radians(units::si::Metres<T>(length().v() / T(2))));
  }

  /// Calculate the geodesic point and pole at the arc length along the
  /// geodesic.
  /// @param gc_length the great circle length on the auxiliary sphere, in
  /// radians.
  /// @return the point and pole projected onto the auxiliary sphere at at
  /// gc_length.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_point_and_pole(const Radians<T> gc_length) const
      -> std::tuple<V, V> {
    auto point{a()};
    auto pole{vector::calculate_pole(beta_, lon_, azi_)};

    if (gc_length.abs().v() < great_circle::MIN_VALUE<T>)
      return {point, pole};

    const Angle<T> length{Angle<T>(gc_length)};
    const Angle<T> beta{arc_beta(length)};
    const Angle<T> lon{arc_longitude(gc_length)};
    point = vector::to_point(beta, lon);

    // if point is on a meridional geodesic, use auxiliary sphere point and pole
    if (azi0_.sin().abs().v() < great_circle::MIN_VALUE<T>)
      return {point, pole};

    // Note: point cannot be at North pole, since it is not on a meridional
    // geodesic. Use Karney's method to calculate azimuth.
    const Angle<T> sigma_sum{sigma1_ + length};
    const auto azimuth{Angle<T>::from_y_x(
        azi0_.sin().v(), azi0_.cos().v() * sigma_sum.cos().v())};
    pole = vector::calculate_pole(beta, lon, azimuth);
    return {point, pole};
  }

  /// Calculate along and across track distances to a position from a geodesic.
  /// @tparam MAX_ITERATIONS the maximum number of iterations, default 10.
  /// @param beta the parametric latitude of the position
  /// @param lon the longitude of the position
  /// @param precision the required precision, in Radians
  ///
  /// @return the along and across track distances to the position in `Radians`.
  template <unsigned MAX_ITERATIONS = 10u>
  [[nodiscard("Pure Function")]]
  constexpr auto calculate_arc_atd_and_xtd(const Angle<T> beta, Angle<T> lon,
                                           Radians<T> precision) const
      -> std::tuple<Radians<T>, Radians<T>, unsigned> {
    // calculate the position as a point on the auxiliary sphere
    const V point{vector::to_point(beta, lon)};
    // calculate the start point and pole on the auxiliary sphere
    const auto [a, pole]{arc_point_and_pole(Radians(T()))};
    const auto gc_d{great_circle::e2gc_distance(vector::distance(a, point))};

    // if the point is close to the start point of the GeodesicSegment
    if (gc_d.v() < precision.v())
      return {Radians<T>(0), Radians<T>(0), 0};

    auto [atd, xtd]{vector::calculate_atd_and_xtd(a, pole, point)};
    auto iterations{1u};
    while (iterations < MAX_ITERATIONS) {
      const Angle<T> atd_angle{atd};
      const auto beta_x{arc_beta(atd_angle)};
      const auto lon_x{arc_longitude(atd)};
      const auto azi_x{arc_azimuth(atd_angle)};

      // calculate the geodesic azimuth and length to the point from the
      // GeodesicSegment position at atd
      const auto [azi_p, length, _azi_end_, p_]{aux_sphere_azimuths_length(
          beta_x, beta, lon - lon_x, Radians<T>(great_circle::MIN_VALUE<T>),
          ellipsoid_)};
      const auto delta_azi{azi_x - azi_p};
      const Radians<T> delta_atd(
          trig::spherical_cosine_rule(delta_azi.cos().v(), length.v()));
      atd = atd + delta_atd;
      xtd = length;

      if (delta_atd.abs().v() < precision.v()) {
        break;
      }

      ++iterations;
    }

    // get the cross track distance (and sign) at the along track distance
    if (xtd.v() < precision.v()) {
      xtd = Radians(T());
    } else {
      const auto [_a, pole_x]{arc_point_and_pole(atd)};
      const auto sign{pole_x.dot(point)};
      xtd = Radians<T>(std::copysign(xtd.v(), sign));
    }

    return {atd, xtd, iterations};
  }

  /// Calculate along and across track distances to a position from a geodesic.
  /// @tparam MAX_ITERATIONS the maximum number of iterations, default 10.
  /// @param position the position as a `LatLong`
  /// @param precision_m the required precision, in Metres
  ///
  /// @return the along and across track distances to the position in `Metres`.
  template <unsigned MAX_ITERATIONS = 10u>
  [[nodiscard("Pure Function")]]
  constexpr auto calculate_atd_and_xtd(const LatLong<T> position,
                                       units::si::Metres<T> precision_m) const
      -> std::tuple<units::si::Metres<T>, units::si::Metres<T>, unsigned> {
    // convert precision to Radians
    const Radians<T> precision{precision_m.v() / ellipsoid_.a().v()};

    // calculate the parametric latitude and longitude of the position
    const Angle<T> beta{
        ellipsoid_.calculate_parametric_latitude(Angle<T>(position.lat()))};
    const Angle<T> lon(position.lon());

    const auto [atd, xtd, iterations]{
        calculate_arc_atd_and_xtd<MAX_ITERATIONS>(beta, lon, precision)};

    // calculate the parametric latitude and azimuth at the abeam point
    const Angle<T> atd_angle{atd};
    const Angle<T> beta_x{arc_beta(atd_angle)};
    const Angle<T> alpha{arc_azimuth(atd_angle).quarter_turn_ccw()};
    return {radians_to_metres(atd),
            convert_radians_to_metres(beta_x, alpha, xtd, ellipsoid_),
            iterations};
  }
};

} // namespace ellipsoid
} // namespace via
