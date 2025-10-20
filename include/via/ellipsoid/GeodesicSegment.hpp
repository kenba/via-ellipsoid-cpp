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
#include <via/angle.hpp>
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
  /// Great circle arc distance to the first Equator crossing.
  Angle<T> sigma1_ = Angle<T>(T(), T(1));
  /// Great circle arc length on the auxiliary sphere in radians.
  Radians<T> arc_length_ = Radians<T>(T());
  /// The half width of a Geodesic Rectangle in metres.
  units::si::Metres<T> half_width_ = units::si::Metres<T>(T());
  /// Integration constant: epsilon, derived from Clairaut's constant.
  T eps_ = T();
  T a3c_ = T(); ///< constant used to convert geodesic/great circle longitudes.
  /// Start parameter for geodesic/great circle distance differences.
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
  /// @param half_width the GeodesicSegment half width in Metres, default zero.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const Angle<T> beta, const Angle<T> lon, const Angle<T> azimuth,
      const Radians<T> arc_length,
      const units::si::Metres<T> half_width = units::si::Metres<T>(T()),
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : beta_{beta}, lon_{lon}, azi_{azimuth},
        // Calculate the azimuth at the first Equator crossing
        azi0_(trig::UnitNegRange(azi_.sin().v() * beta_.cos().v()),
              trig::swap_sin_cos(
                  trig::UnitNegRange(azi_.sin().v() * beta_.cos().v()))),
        // Calculate the distance to the first Equator crossing
        sigma1_{Angle<T>::from_y_x(beta_.sin().v(),
                                   beta.cos().v() * azi_.cos().v())},
        arc_length_{arc_length}, half_width_{half_width},
        // Calculate eps for calculating coefficients
        eps_(ellipsoid.calculate_epsilon(azi0_.sin())),
        a3c_(ellipsoid.calculate_a3c(azi0_.sin(), eps_)),
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
  /// @param half_width the GeodesicSegment half width in Metres, default zero.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const LatLong<T> &a, const Angle<T> &azimuth,
      const Radians<T> arc_length = Radians(T()),
      const units::si::Metres<T> half_width = units::si::Metres<T>(T()),
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : GeodesicSegment(
            ellipsoid.calculate_parametric_latitude(Angle<T>(a.lat())),
            Angle<T>(a.lon()), azimuth, arc_length, half_width, ellipsoid) {}

  /// Construct a GeodesicSegment from a start point, azimuth and length in
  /// Metres.
  /// @param a the geodetic start point.
  /// @param azimuth the azimuth at the start point.
  /// @param length the length in Metres.
  /// @param half_width the GeodesicSegment half width in Metres, default zero.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const LatLong<T> &a, const Angle<T> &azimuth,
      const units::si::Metres<T> length,
      const units::si::Metres<T> half_width = units::si::Metres<T>(T()),
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : GeodesicSegment(
            ellipsoid.calculate_parametric_latitude(Angle<T>(a.lat())),
            Angle<T>(a.lon()), azimuth, Radians(T()), half_width, ellipsoid) {
    arc_length_ = metres_to_radians(length);
  }

  /// Construct a GeodesicSegment from a start point and azimuth, arc length
  /// tuple.
  /// @param a the geodetic point.
  /// @param azimuth_length the azimuth, arc length tuple.
  /// @param half_width the GeodesicSegment half width in Metres, default zero.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const LatLong<T> &a,
      const std::tuple<Angle<T>, Radians<T>, Angle<T>, unsigned> azimuth_length,
      const units::si::Metres<T> half_width = units::si::Metres<T>(T()),
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : GeodesicSegment(a, std::get<0>(azimuth_length),
                        std::get<1>(azimuth_length), half_width, ellipsoid) {}

  /// Construct a GeodesicSegment between a pair of points.
  ///
  /// @pre a.is_valid() && b.s_valid()
  /// @pre tolerance >= epsilon
  ///
  /// @param a, b the start and finish points in geodetic coordinates.
  /// @param half_width the GeodesicSegment half width in Metres, default zero.
  /// @param tolerance the tolerance to perform the calculation to in Radians,
  /// default great_circle::MIN_VALUE.
  /// @param ellipsoid a const reference to the underlying Ellipsoid, default
  /// wgs84.
  constexpr GeodesicSegment(
      const LatLong<T> &a, const LatLong<T> &b,
      const units::si::Metres<T> half_width = units::si::Metres<T>(T()),
      const Radians<T> tolerance = Radians<T>(great_circle::MIN_VALUE<T>),
      const Ellipsoid<T> &ellipsoid = Ellipsoid<T>::wgs84())
      : GeodesicSegment(
            a, calculate_azimuths_arc_length(a, b, tolerance, ellipsoid),
            half_width, ellipsoid) {}

  /// Test whether a `GeodesicSegment` is valid:
  /// @invariant 0° <= latitude <= 90°
  /// @invariant 0 <= arc_length <= π
  [[nodiscard("Pure Function")]]
  constexpr auto is_valid() const noexcept -> bool {
    return beta_.cos().v() >= T() && T() <= arc_length_.v() &&
           arc_length_.v() <= trig::PI<T> && half_width_.v() >= T();
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

  /// Set the `arc_length` of a `GeodesicSegment`
  /// @param arc_length the great circle arc length of the `GeodesicSegment`.
  constexpr auto set_arc_length(units::si::Metres<T> arc_length) noexcept
      -> void {
    arc_length_ = arc_length;
  }

  /// Accessor for the arc length on the auxiliary sphere in radians.
  /// @return The arc length on the auxiliary sphere in radians.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_length() const noexcept -> Radians<T> {
    return arc_length_;
  }

  /// Set the `half_width` of a `GeodesicSegment`
  /// @param half_width of the `GeodesicSegment`.
  constexpr auto set_half_width(units::si::Metres<T> half_width) noexcept
      -> void {
    half_width_ = half_width;
  }

  /// Accessor for the half width in metres.
  [[nodiscard("Pure Function")]]
  constexpr auto half_width() const noexcept -> units::si::Metres<T> {
    return half_width_;
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
  /// @param distance the distance along the GeodesicSegment in metres.
  /// @return the distance along the great circle arc in radians.
  [[nodiscard("Pure Function")]]
  constexpr auto metres_to_radians(const units::si::Metres<T> distance) const
      -> Radians<T> {
    if (std::abs(distance.v()) < great_circle::MIN_VALUE<T>)
      return Radians(T());

    const auto a1((evaluate_A1<T>(eps_) + T(1)));
    const auto tau12{Radians<T>(distance.v() / (ellipsoid_.b().v() * a1))};
    const auto tau_sum{Angle<T>(b11_ + tau12)};
    const auto c1p{evaluate_coeffs_C1p<T>(eps_)};
    const auto b12{sin_cos_series(sigma1_ + Angle<T>(tau_sum), c1p)};

    return tau12 + b12 + b11_;
  }

  /// Convert a distance in radians on the auxiliary sphere to metres
  /// on the ellipsoid.
  /// @param arc_distance the distance along the great circle arc in Radians.
  /// @param sigma the arc_distance as an Angle.
  /// @return the distance along the GeodesicSegment in metres.
  [[nodiscard("Pure Function")]]
  constexpr auto radians_to_metres(const Radians<T> arc_distance,
                                   const Angle<T> sigma) const
      -> units::si::Metres<T> {
    // Calculate Great circle distance from Northward Equator crossing.
    const Angle<T> sigma_sum{sigma1_ + sigma};
    const auto c1{evaluate_coeffs_C1<T>(eps_)};
    const auto b12{sin_cos_series(sigma_sum, c1)};
    const auto a1((evaluate_A1<T>(eps_) + T(1)));

    return units::si::Metres<T>(ellipsoid().b().v() * a1 *
                                (arc_distance + b12 - b11_).v());
  }

  /// Accessor for the length of the GeodesicSegment in metres.
  /// @return The length of the GeodesicSegment in metres.
  [[nodiscard("Pure Function")]]
  constexpr auto length() const -> units::si::Metres<T> {
    return radians_to_metres(arc_length_, Angle<T>(arc_length_));
  }

  /// Calculate the parametric latitude at the great circle distance.
  /// @param sigma the arc distance on the auxiliary sphere as an Angle.
  /// @return the parametric latitude at sigma.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_beta(const Angle<T> sigma) const -> Angle<T> {
    return great_circle::calculate_latitude(beta_, azi_, sigma);
  }

  /// Calculate the geodetic latitude at the great circle arc distance.
  /// @param arc_distance the great circle arc distance on the auxiliary
  /// sphere.
  /// @return the geodetic latitude of the position at arc_distance.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_latitude(const Radians<T> arc_distance) const -> Angle<T> {
    return ellipsoid_.calculate_geodetic_latitude(
        arc_beta(Angle<T>(arc_distance)));
  }

  /// Calculate the geodetic latitude at distance along the geodesic.
  /// @param distance the distance along the geodesic segment, in metres.
  /// @return the geodetic latitude at distance.
  [[nodiscard("Pure Function")]]
  constexpr auto latitude(const units::si::Metres<T> distance) const
      -> Angle<T> {
    return arc_latitude(metres_to_radians(distance));
  }

  /// Calculate the azimuth along the arc at great circle distance sigma.
  /// @param sigma the arc distance on the auxiliary sphere as an Angle.
  /// @return the azimuth at `sigma`.
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

  /// Calculate the azimuth at the distance along the geodesic.
  /// @param distance the distance along the geodesic segment, in metres.
  /// @return the azimuth of the geodesic segment at distance.
  [[nodiscard("Pure Function")]]
  constexpr auto
  azimuth(const units::si::Metres<T> distance = units::si::Metres(T())) const
      -> Angle<T> {
    const Angle<T> sigma{metres_to_radians(distance)};
    return arc_azimuth(sigma);
  }

  /// Calculate the geodesic longitude difference at arc distance
  /// along the auxiliary sphere.
  /// @param arc_distance the great circle arc distance on the auxiliary
  /// sphere.
  /// @param sigma the arc_distance as an Angle.
  /// @return the longitude difference from the start point.
  [[nodiscard("Pure Function")]]
  constexpr auto delta_longitude(const Radians<T> arc_distance,
                                 const Angle<T> sigma) const -> Angle<T> {
    if (arc_distance.abs().v() < great_circle::MIN_VALUE<T>)
      return Angle<T>();

    // The great circle distance from Northward Equator crossing.
    const Angle<T> sigma_sum(sigma1_ + sigma);

    // The longitude difference on the auxiliary sphere.
    const auto omega1{Angle<T>::from_y_x(azi0_.sin().v() * beta_.sin().v(),
                                         beta_.cos().v() * azi_.cos().v())};
    const auto omega2{Angle<T>::from_y_x(azi0_.sin().v() * sigma_sum.sin().v(),
                                         sigma_sum.cos().v())};
    const Angle<T> omega12{omega2 - omega1};

    const auto c3{ellipsoid_.calculate_c3y(eps_)};
    const auto b31(sin_cos_series(sigma1_, c3));
    const auto b32(sin_cos_series(sigma_sum, c3));

    return omega12 -
           Angle<T>(Radians<T>{a3c_ * (arc_distance + (b32 - b31)).v()});
  }

  /// Calculate the geodesic longitude at the great circle length along
  /// the auxiliary sphere.
  /// @param arc_distance the great circle arc distance on the auxiliary
  /// sphere.
  /// @return the geodesic longitude at arc_distance.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_longitude(const Radians<T> arc_distance) const
      -> Angle<T> {
    const Angle<T> sigma{arc_distance};
    return lon_ + delta_longitude(arc_distance, sigma);
  }

  /// Calculate the geodesic longitude at the length along the geodesic.
  /// @param distance the distance along the geodesic segment, in metres.
  /// @return the geodetic longitude as an Angle.
  [[nodiscard("Pure Function")]]
  constexpr auto longitude(const units::si::Metres<T> distance) const
      -> Angle<T> {
    return arc_longitude(metres_to_radians(distance));
  }

  /// Calculate the parametric latitude and azimuth at the arc distance.
  /// @param arc_distance the arc distance on the auxiliary sphere in Radians.
  /// @return the parametric latitude and azimuth at arc distance.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_beta_long(const Radians<T> arc_distance) const
      -> std::tuple<Angle<T>, Angle<T>> {
    const Angle<T> sigma{arc_distance};
    const Angle<T> beta{arc_beta(sigma)};
    const Angle<T> lon{lon_ + delta_longitude(arc_distance, sigma)};

    return {beta, lon};
  }

  /// Calculate the geodesic `LatLong` at the great circle length along
  /// the auxiliary sphere.
  /// @param arc_distance the great circle arc distance on the auxiliary
  /// sphere.
  /// @return the `LatLong` of the geodesic position at arc_distance.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_lat_long(const Radians<T> arc_distance) const
      -> LatLong<T> {
    const auto [beta, lon]{arc_beta_long(arc_distance)};
    return LatLong<T>(ellipsoid_.calculate_geodetic_latitude(beta).to_degrees(),
                      lon.to_degrees());
  }

  /// Calculate the geodesic `LatLong` at a distance along the
  /// `GeodesicSegment`.
  /// @param distance the distance along the geodesic segment, in metres.
  /// @return the `LatLong` of the geodesic position at distance.
  [[nodiscard("Pure Function")]]
  constexpr auto lat_long(const units::si::Metres<T> distance) const
      -> LatLong<T> {
    return arc_lat_long(metres_to_radians(distance));
  }

  /// Calculate the parametric latitude, longitude and azimuth at the arc
  /// distance.
  /// @param arc_distance the arc distance on the auxiliary sphere in Radians.
  /// @return the parametric latitude, longitude and azimuth at arc distance.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_angles(const Radians<T> arc_distance) const
      -> std::tuple<Angle<T>, Angle<T>, Angle<T>> {
    const Angle<T> sigma{arc_distance};
    const Angle<T> beta{arc_beta(sigma)};
    const Angle<T> lon{lon_ + delta_longitude(arc_distance, sigma)};
    const Angle<T> azimuth{arc_azimuth(sigma)};

    return {beta, lon, azimuth};
  }

  /// Calculate the point on the auxiliary sphere at the
  /// great circle arc length.
  /// @param arc_distance the great circle arc distance on the auxiliary
  /// sphere.
  /// @return the point on the auxiliary sphere at arc_distance.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_point(const Radians<T> arc_distance) const -> V {
    if (arc_distance.abs().v() < great_circle::MIN_VALUE<T>)
      return a();

    const auto [beta, lon]{arc_beta_long(arc_distance)};
    return vector::to_point(beta, lon);
  }

  /// Calculate the point on the auxiliary sphere at the mid point of the
  /// `GeodesicSegment`.
  /// @return mid point
  [[nodiscard("Pure Function")]]
  constexpr auto mid_point() const -> V {
    return arc_point(metres_to_radians(length().half()));
  }

  /// Calculate the geodesic pole at the great circle arc distance.
  /// @param arc_distance the great circle distance on the auxiliary sphere
  /// @return the pole projected onto the auxiliary sphere at at arc_distance.
  constexpr auto arc_pole(const Radians<T> arc_distance) const -> V {
    // if point is on a meridional geodesic, use auxiliary sphere point and
    // pole
    if (azi0_.sin().abs().v() < great_circle::MIN_VALUE<T>) {
      return vector::calculate_pole(beta_, lon_, azi_);
    }

    const auto [beta, lon, azimuth]{arc_angles(arc_distance)};
    return vector::calculate_pole(beta, lon, azimuth);
  }

  /// Calculate the geodesic point and pole at the arc distance along the
  /// geodesic.
  /// @param arc_distance the great circle distance on the auxiliary sphere
  /// @return the point and pole projected onto the auxiliary sphere at at
  /// arc_distance.
  [[nodiscard("Pure Function")]]
  constexpr auto arc_point_and_pole(const Radians<T> arc_distance) const
      -> std::tuple<V, V> {
    const auto [beta, lon, azimuth]{arc_angles(arc_distance)};
    const auto pole{(azi0_.sin().abs().v() < great_circle::MIN_VALUE<T>)
                        ? vector::calculate_pole(beta_, lon_, azi_)
                        : vector::calculate_pole(beta, lon, azimuth)};
    return {vector::to_point(beta, lon), pole};
  }

  /// The reverse `GeodesicSegment` from end to start.
  ///
  /// @return the reverse `GeodesicSegment` from end to start.
  [[nodiscard("Pure Function")]]
  constexpr auto reverse() const -> GeodesicSegment<T> {
    const Angle<T> sigma(arc_length_);
    return GeodesicSegment<T>(
        arc_beta(sigma), lon_ + delta_longitude(arc_length_, sigma),
        arc_azimuth(sigma).opposite(), arc_length_, half_width_, ellipsoid_);
  }

  /// Calculate along and across track distances to a position from a
  /// geodesic.
  /// @tparam MAX_ITERATIONS the maximum number of iterations, default 10.
  /// @param beta the parametric latitude of the position
  /// @param lon the longitude of the position
  /// @param precision the required precision, in Radians
  ///
  /// @return the along and across track distances to the position in
  /// `Radians`.
  template <unsigned MAX_ITERATIONS = 10u>
  [[nodiscard("Pure Function")]]
  constexpr auto calculate_arc_atd_and_xtd(const Angle<T> beta, Angle<T> lon,
                                           const Radians<T> precision) const
      -> std::tuple<Radians<T>, Radians<T>, unsigned> {
    // calculate the position as a point on the auxiliary sphere
    const V point{vector::to_point(beta, lon)};

    // if the point is close to the start point of the GeodesicSegment
    const auto gc_d{great_circle::e2gc_distance(vector::distance(a(), point))};
    if (gc_d.v() < precision.v())
      return {Radians<T>(0), Radians<T>(0), 0};

    const auto pole{vector::calculate_pole(beta_, lon_, azi_)};
    auto [atd, xtd]{vector::calculate_atd_and_xtd(a(), pole, point)};
    auto iterations{1u};
    while (iterations < MAX_ITERATIONS) {
      const auto [beta_x, lon_x, azi_x]{arc_angles(atd)};

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
      const auto pole_x{arc_pole(atd)};
      const auto sign{pole_x.dot(point)};
      xtd = Radians<T>(std::copysign(xtd.v(), sign));
    }

    return {atd, xtd, iterations};
  }

  /// Calculate the shortest geodesic distance of point from the
  /// `GeodesicSegment`.
  /// @tparam MAX_ITERATIONS the maximum number of iterations, default 10.
  /// @param beta the parametric latitude of the point.
  /// @param lon the longitude of the point.
  /// @param precision the required precision, in `Radians`.
  ///
  /// @return shortest distance of the point from the `GeodesicSegment` in
  /// `Metres`.
  template <unsigned MAX_ITERATIONS = 10u>
  [[nodiscard("Pure Function")]]
  constexpr auto
  calculate_sphere_shortest_distance(const Angle<T> beta, Angle<T> lon,
                                     const Radians<T> precision) const
      -> units::si::Metres<T> {
    const auto [atd, xtd, _]{
        calculate_arc_atd_and_xtd<MAX_ITERATIONS>(beta, lon, precision)};

    // if the position is beside the geodesic segment
    if (vector::intersection::is_alongside(atd, arc_length_, precision)) {
      if (xtd.abs().v() < precision.v()) {
        return units::si::Metres<T>(0);
      } else {
        // convert cross track distance to Metres
        const auto [beta_x, lon, azi]{arc_angles(atd)};
        const Angle<T> alpha{azi.quarter_turn_ccw()};
        const auto distance{
            convert_radians_to_metres(beta_x, alpha, xtd, ellipsoid_)};
        // return the abs cross track distance in Metres
        return units::si::Metres<T>(std::abs(distance.v()));
      }
    } else {
      // adjust atd to measure the distance from the centre of the Arc
      const auto atd_centre{atd - arc_length_.half()};
      if (std::signbit(atd_centre.v())) {
        // calculate the geodesic distance from the start of the segment
        const Angle<T> delta_long{lon - lon_};
        const auto [alpha, distance, _1, _2] = aux_sphere_azimuths_length(
            beta_, beta, delta_long, precision, ellipsoid_);
        return convert_radians_to_metres(beta_, alpha, distance, ellipsoid_);
      } else {
        // calculate the geodesic distance from the end of the segment
        const auto [beta_x, lon_x, _azi]{arc_angles(arc_length_)};
        const Angle<T> delta_long{lon - lon_x};

        const auto [alpha, distance, _1, _2] = aux_sphere_azimuths_length(
            beta_x, beta, delta_long, precision, ellipsoid_);
        return convert_radians_to_metres(beta_x, alpha, distance, ellipsoid_);
      }
    }
  }

  /// Calculate along and across track distances to a position from a
  /// geodesic.
  /// @tparam MAX_ITERATIONS the maximum number of iterations, default 10.
  /// @param position the position as a `LatLong`
  /// @param precision_m the required precision, in Metres
  ///
  /// @return the along and across track distances to the position in
  /// `Metres`.
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
    return {radians_to_metres(atd, atd_angle),
            convert_radians_to_metres(beta_x, alpha, xtd, ellipsoid_),
            iterations};
  }

  /// Calculate the shortest geodesic distance of point from the
  /// `GeodesicSegment`.
  /// @tparam MAX_ITERATIONS the maximum number of iterations, default 10.
  /// @param position the position as a `LatLong`
  /// @param precision_m the required precision, in `Metres`.
  ///
  /// @return shortest distance of the point from the `GeodesicSegment` in
  /// `Metres`.
  template <unsigned MAX_ITERATIONS = 10u>
  [[nodiscard("Pure Function")]]
  constexpr auto shortest_distance(const LatLong<T> position,
                                   units::si::Metres<T> precision_m) const
      -> units::si::Metres<T> {
    // calculate the parametric latitude and longitude of the position
    const Angle<T> beta{
        ellipsoid_.calculate_parametric_latitude(Angle<T>(position.lat()))};
    const Angle<T> lon(position.lon());
    // convert precision to Radians
    const Radians<T> precision{precision_m.v() / ellipsoid_.a().v()};
    return calculate_sphere_shortest_distance<MAX_ITERATIONS>(beta, lon,
                                                              precision);
  }
};

} // namespace ellipsoid
} // namespace via
