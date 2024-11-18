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
/// @file coefficients.hpp
/// @brief Contains the via::ellipsoid coefficients and series expansion
/// functions.
//////////////////////////////////////////////////////////////////////////////
#include <array>
#include <via/angle.hpp>

namespace via {
namespace ellipsoid {
/// The scale factor A1.
/// @see CFF Karney, Geodesics on an ellipsoid of revolution: Eq. 48,
/// https://arxiv.org/pdf/1102.1215.pdf.
/// @tparam T a floating point type, e.g.: double or long double.
/// @param eps epsilon the integration variable derived from Clairaut's
/// constant.
/// @return A1
template <typename T,
          size_t order = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto evaluate_A1(const T eps) noexcept -> T {
  const T eps2{eps * eps};

  if constexpr (std::is_same<T, long double>::value) {
    const T t{eps2 * (eps2 * (eps2 * (25 * eps2 + 64) + 256) + 4096) /
              T(16384)};
    return (t + eps) / (T(1) - eps);
  } else {
    const T t{eps2 * (eps2 * (eps2 + 4) + 64) / T(256)};
    return (t + eps) / (T(1) - eps);
  }
}

/// The scale factor A2.
/// @see CFF Karney, Geodesics on an ellipsoid of revolution: Eq. 50,
/// https://arxiv.org/pdf/1102.1215.pdf.
/// @tparam T a floating point type, e.g.: double or long double.
/// @param eps epsilon the integration variable derived from Clairaut's
/// constant.
/// @return A2
template <typename T,
          size_t order = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto evaluate_A2(const T eps) noexcept -> T {
  const T eps2{eps * eps};
  if constexpr (std::is_same<T, long double>::value) {
    const T t{eps2 * (eps2 * ((-375 * eps2 - 704) * eps2 - 1792) - 12288) /
              T(16384)};
    return (t - eps) / (1 + eps);
  } else {
    const T t{eps2 * ((-11 * eps2 - 28) * eps2 - 192) / T(256)};
    return (t - eps) / (1 + eps);
  }
}

/// The scale factor A3.
/// @see CFF Karney, Geodesics on an ellipsoid of revolution: Eq. 52,
/// https://arxiv.org/pdf/1102.1215.pdf.
/// @tparam T a floating point type, e.g.: double or long double.
/// @param n the third flattening of the ellipsoid.
/// @return A3
template <typename T,
          size_t order = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto evaluate_coeffs_A3(const T n) noexcept -> std::array<T, order> {
  if constexpr (std::is_same<T, long double>::value)
    return std::array<T, order>{T(1),
                                (n - 1) / T(2),
                                (n * (3 * n - 1) - 2) / T(8),
                                (n * (n * (5 * n - 1) - 3) - 1) / T(16),
                                (n * ((-5 * n - 20) * n - 4) - 6) / T(128),
                                ((-5 * n - 10) * n - 6) / T(256),
                                (-15 * n - 20) / T(1024),
                                -25 / T(2048)};
  else
    return std::array<T, order>{T(1),
                                (n - 1) / T(2),
                                (n * (3 * n - 1) - 2) / T(8),
                                ((-n - 3) * n - 1) / T(16),
                                (-2 * n - 3) / T(64),
                                -3 / T(128)};
}

/// The coefficients C1[l] in the Fourier expansion of B1.
/// @see CFF Karney, Geodesics on an ellipsoid of revolution: Eq. 49,
/// https://arxiv.org/pdf/1102.1215.pdf.
/// @tparam T a floating point type, e.g.: double or long double.
/// @param eps epsilon the integration variable derived from Clairaut's
/// constant.
/// @return C1 coefficients array.
template <typename T,
          size_t order = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto evaluate_coeffs_C1(const T eps) noexcept
    -> std::array<T, order + 1> {
  const T eps2{eps * eps};
  const T eps4{(eps2 * eps) * eps}; // Note: not the same as eps2 * eps2!
  const T eps6{(eps4 * eps) * eps};

  if constexpr (std::is_same<T, long double>::value)
    return std::array<T, order + 1>{
        T(0),
        eps * (eps2 * (eps2 * (19 * eps2 - 64) + 384) - 1024) / T(2048),
        eps2 * (eps2 * (eps2 * (7 * eps2 - 18) + 128) - 256) / T(4096),
        eps * eps2 * ((72 - 9 * eps2) * eps2 - 128) / T(6144),
        eps4 * ((96 - 11 * eps2) * eps2 - 160) / T(16384),
        eps * eps4 * (35 * eps2 - 56) / T(10240),
        eps6 * (9 * eps2 - 14) / T(4096),
        eps * eps6 * -33 / T(14336),
        eps * (eps * eps6) * -429 / T(262144)};
  else
    return std::array<T, order + 1>{T(0),
                                    eps * ((6 - eps2) * eps2 - 16) / T(32),
                                    eps2 * ((64 - 9 * eps2) * eps2 - 128) /
                                        T(2048),
                                    eps * eps2 * (9 * eps2 - 16) / T(768),
                                    eps4 * (3 * eps2 - 5) / T(512),
                                    eps * eps4 * (-7 / T(1280)),
                                    eps6 * (-7 / T(2048))};
}

/// The coefficients C1p[l] in the Fourier expansion of B1p.
/// @see CFF Karney, Geodesics on an ellipsoid of revolution: Eq. 58,
/// https://arxiv.org/pdf/1102.1215.pdf.
/// @tparam T a floating point type, e.g.: double or long double.
/// @param eps epsilon the integration variable derived from Clairaut's
/// constant.
/// @return C1p coefficients array.
template <typename T,
          size_t order = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto evaluate_coeffs_C1p(const T eps) noexcept
    -> std::array<T, order> {
  const T eps2{eps * eps};
  const T eps4{(eps2 * eps) * eps}; // Note: not the same as eps2 * eps2!

  if constexpr (std::is_same<T, long double>::value) {
    const T eps6{(eps4 * eps) * eps};
    return std::array<T, order>{
        T(0),
        eps * (eps2 * ((9840 - 4879 * eps2) * eps2 - 20736) + 36864) / T(73728),
        eps2 * (eps2 * (4005 * eps2 - 4736) + 3840) / T(12288),
        eps * eps2 * (eps2 * (8703 * eps2 - 7200) + 3712) / T(12288),
        eps4 * (2695 - 7173 * eps2) / T(7680),
        (eps * eps4) * (41604 - 141115 * eps2) / T(92160),
        eps6 * 38081 / T(61440),
        (eps * eps6) * 459485 / T(516096)};
  } else
    return std::array<T, order>{T(0),
                                eps * (eps2 * (205 * eps2 - 432) + 768) /
                                    T(1536),
                                eps2 * (30 - 37 * eps2) / T(96),
                                eps * eps2 * (116 - 225 * eps2) / T(384),
                                eps4 * 539 / T(1536),
                                (eps * eps4) * 3467 / T(7680)};
}

/// The coefficients C2[l] in the Fourier expansion of B2.
template <typename T,
          size_t order = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto evaluate_coeffs_C2(const T eps) noexcept
    -> std::array<T, order + 1> {
  const T eps2{eps * eps};
  const T eps4{(eps2 * eps) * eps}; // Note: not the same as eps2 * eps2!
  const T eps6{(eps4 * eps) * eps};

  if constexpr (std::is_same<T, long double>::value)
    return std::array<T, order + 1>{
        T(0),
        eps * (eps2 * (eps2 * (41 * eps2 + 64) + 128) + 1024) / T(2048),
        eps2 * (eps2 * (eps2 * (47 * eps2 + 70) + 128) + 768) / T(4096),
        eps * eps2 * (eps2 * (69 * eps2 + 120) + 640) / T(6144),
        eps4 * (eps2 * (133 * eps2 + 224) + 1120) / T(16384),
        eps * eps4 * (105 * eps2 + 504) / T(10240),
        eps6 * (33 * eps2 + 154) / T(4096),
        eps * eps6 * 429 / T(14336),
        eps * (eps * eps6) * 6435 / T(262144)};
  else
    return std::array<T, order + 1>{T(0),
                                    eps * (eps2 * (eps2 + 2) + 16) / T(32),
                                    eps2 * (eps2 * (35 * eps2 + 64) + 384) /
                                        T(2048),
                                    eps * eps2 * (15 * eps2 + 80) / T(768),
                                    eps4 * (7 * eps2 + 35) / T(512),
                                    eps * eps4 * 63 / T(1280),
                                    eps6 * 77 / T(2048)};
}

/// The coefficients C3x[l] in the Fourier expansion of C3.
/// @see CFF Karney, Geodesics on an ellipsoid of revolution: Eq. 53,
/// https://arxiv.org/pdf/1102.1215.pdf.
/// @tparam T a floating point type, e.g.: double or long double.
/// @param n the third flattening of the ellipsoid.
/// @return C3x coefficients array.
template <typename T,
          size_t order = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto evaluate_coeffs_C3x(const T n) noexcept
    -> std::array<T, order *(order - 1) / 2> {
  const T n2(n * n);

  if constexpr (std::is_same<T, long double>::value)
    return std::array<T, order *(order - 1) / 2>{
        (1 - n) / T(4),
        (1 - n2) / T(8),
        (n * ((-5 * n - 1) * n + 3) + 3) / T(64),
        (n * ((2 - 2 * n) * n + 2) + 5) / T(128),
        (n * (3 * n + 11) + 12) / T(512),
        (10 * n + 21) / T(1024),
        243 / T(16384),
        ((n - 3) * n + 2) / T(32),
        (n * (n * (2 * n - 3) - 2) + 3) / T(64),
        (n * ((-6 * n - 9) * n + 2) + 6) / T(256),
        ((1 - 2 * n) * n + 5) / T(256),
        (69 * n + 108) / T(8192),
        187 / T(16384),
        (n * ((5 - n) * n - 9) + 5) / T(192),
        (n * (n * (10 * n - 6) - 10) + 9) / T(384),
        ((-77 * n - 8) * n + 42) / T(3072),
        (12 - n) / T(1024),
        139 / T(16384),
        (n * ((20 - 7 * n) * n - 28) + 14) / T(1024),
        ((-7 * n - 40) * n + 28) / T(2048),
        (72 - 43 * n) / T(8192),
        127 / T(16384),
        (n * (75 * n - 90) + 42) / T(5120),
        (9 - 15 * n) / T(1024),
        99 / T(16384),
        (44 - 99 * n) / T(8192),
        99 / T(16384),
        429 / T(114688)};
  else
    return std::array<T, order *(order - 1) / 2>{
        (1 - n) / T(4),
        (1 - n2) / T(8),
        (n * ((-5 * n - 1) * n + 3) + 3) / T(64),
        (n * ((2 - 2 * n) * n + 2) + 5) / T(128),
        (n * (3 * n + 11) + 12) / T(512),
        ((n - 3) * n + 2) / T(32),
        (n * (n * (2 * n - 3) - 2) + 3) / T(64),
        (n * ((-6 * n - 9) * n + 2) + 6) / T(256),
        ((1 - 2 * n) * n + 5) / T(256),
        (n * ((5 - n) * n - 9) + 5) / T(192),
        (n * (n * (10 * n - 6) - 10) + 9) / T(384),
        ((-77 * n - 8) * n + 42) / T(3072),
        (n * ((20 - 7 * n) * n - 28) + 14) / T(1024),
        ((-7 * n - 40) * n + 28) / T(2048),
        (n * (75 * n - 90) + 42) / T(5120)};
}

/// Evaluate the polynomial in x using Horner's method.
/// @see https://en.wikipedia.org/wiki/Horner%27s_method
/// @pre first < last
/// @tparam T a floating point type, e.g.: double or long double.
/// @param x the variable.
/// @param first, last iterators to the polynomial coefficients.
/// @return the result of evaluating the polynomial.
template <typename T, typename IteratorType>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto polyval(const T x, const IteratorType first,
                       const IteratorType last) -> T {
  Expects(1 <= std::distance(first, last));

  IteratorType it(last);
  auto result(static_cast<T>(*(--it)));
  while (it != first)
    result = x * result + static_cast<T>(*(--it));

  return result;
}

/// Evaluate the polynomial in x using an unrolled version of Horner's method.
/// @pre n > 0
/// @param x the variable.
/// @param coeffs a pointer to the first polynomial coefficient.
/// @param n the number of polynomial coefficients to evaluate.
/// @return the result of evaluating the polynomial.
template <typename T, typename S>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto polyval_n(const T x, const S *coeffs, size_t n) -> T {
  Expects(0 < n);

  T result(0);
  switch (n) {
  case 8u:
    result = static_cast<T>(*(coeffs + 7)) + x * result;
    [[fallthrough]];
  case 7u:
    result = static_cast<T>(*(coeffs + 6)) + x * result;
    [[fallthrough]];
  case 6u:
    result = static_cast<T>(*(coeffs + 5)) + x * result;
    [[fallthrough]];
  case 5u:
    result = static_cast<T>(*(coeffs + 4)) + x * result;
    [[fallthrough]];
  case 4u:
    result = static_cast<T>(*(coeffs + 3)) + x * result;
    [[fallthrough]];
  case 3u:
    result = static_cast<T>(*(coeffs + 2)) + x * result;
    [[fallthrough]];
  case 2u:
    result = static_cast<T>(*(coeffs + 1)) + x * result;
    [[fallthrough]];
  case 1u:
    result = static_cast<T>(*coeffs) + x * result;
    break;

  default:
    result = polyval(x, coeffs, coeffs + n);
  }

  return result;
}

/// Evaluate the polynomial in x.
/// @pre 1 <= coeffs.size()
/// @tparam T a floating point type, e.g.: double or long double.
/// @param x the variable.
/// @param coeffs the polynomial coefficients.
/// @return the result of evaluating the polynomial.
template <typename T, typename Coeffs>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto evaluate_poynomial(const T x, Coeffs const &coeffs) -> T {
  return polyval_n(x, coeffs.data(), coeffs.size());
}

/// The coefficients C3[l] in the Fourier expansion of C3.
/// @see CFF Karney, Geodesics on an ellipsoid of revolution: Eq. 53,
/// https://arxiv.org/pdf/1102.1215.pdf.
/// @tparam T a floating point type, e.g.: double or long double.
/// @param eps epsilon the integration variable derived from Clairaut's
/// constant.
/// @param coeffs2 the polynomial coefficients from evaluate_coeffs_C3x.
/// @return C3 coefficients array.
template <typename T,
          size_t order = std::is_same<T, long double>::value ? 8u : 6u>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto
evaluate_coeffs_C3y(std::array<T, order *(order - 1) / 2> const &coeffs2, T eps)
    -> std::array<T, order> {
  std::array<T, order> coeffs1{0};

  T mult(1);
  int offset(0);
  int index(1);
  // m is order of polynomial in eps.
  for (int m(order - 1); 0 < m; --m, ++index) {
    mult *= eps;
    const auto coeffs2_start(coeffs2.data() + offset);
    coeffs1[index] = mult * polyval_n(eps, coeffs2_start, m);
    offset += m;
  }

  return coeffs1;
}

/// Evaluate the following:
///   y = sum(c[i] * sin(2*i * angle), i, 1, n)
/// using Clenshaw summation.
/// @see https://en.wikipedia.org/wiki/Clenshaw_algorithm
/// @pre 2 <= coeffs.size()
/// @param angle the angle.
/// @param coeffs the polynomial coefficients.
/// @return the result of evaluating the polynomial.
template <typename T, typename Coeffs>
  requires std::floating_point<T>
[[nodiscard("Pure Function")]]
constexpr auto sin_cos_series(const Angle<T> angle, Coeffs const &coeffs)
    -> Radians<T> {
  const Angle<T> angle2x{angle.x2()};
  if (angle2x.sin().abs().v() < std::numeric_limits<T>::epsilon())
    return Radians(T());

  // the Clenshaw ak(theta) parameter, beta(k) = -1
  const T ar{2 * angle2x.cos().v()};

  // Point to last element.
  auto index(coeffs.size() - 1);

  // k1 is the previous element, k0 is the current element
  // If index is odd, k1 = 0, else last coeff
  T k1{(index & 1) ? T() : T(coeffs[index--])};
  T k0{T(coeffs[index--]) + ar * k1};

  // Unroll loop x 2, so accumulators return to their original role.
  while (8u < index) {
    k1 = static_cast<T>(coeffs[index--]) + (ar * k0 - k1);
    k0 = static_cast<T>(coeffs[index--]) + (ar * k1 - k0);
  }

  // Unroll loop
  switch (index) {
  case 8u:
    k1 = static_cast<T>(coeffs[8]) + (ar * k0 - k1);
    k0 = static_cast<T>(coeffs[7]) + (ar * k1 - k0);
    [[fallthrough]];
  case 6u:
    k1 = static_cast<T>(coeffs[6]) + (ar * k0 - k1);
    k0 = static_cast<T>(coeffs[5]) + (ar * k1 - k0);
    [[fallthrough]];
  case 4u:
    k1 = static_cast<T>(coeffs[4]) + (ar * k0 - k1);
    k0 = static_cast<T>(coeffs[3]) + (ar * k1 - k0);
    [[fallthrough]];
  case 2u:
    k1 = static_cast<T>(coeffs[2]) + (ar * k0 - k1);
    k0 = static_cast<T>(coeffs[1]) + (ar * k1 - k0);
  }

  return Radians(angle2x.sin().v() * k0);
}

} // namespace ellipsoid
} // namespace via
