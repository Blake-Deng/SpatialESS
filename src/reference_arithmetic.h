#ifndef SPATIALESS_REFERENCE_ARITHMETIC_HPP
#define SPATIALESS_REFERENCE_ARITHMETIC_HPP

#include <cmath>
#include <array>
#include <cstddef>
#include <limits>
#include <vector>

namespace spatialess {

// Match R's real_mean() accumulation and residual correction, rather than a
// double sum/n. The reference CellChat geometricMean() calls mean(log(x)).
// Adapted from R 4.5 src/main/summary.c (R Core Team, GPL-2-or-later):
// https://github.com/wch/r-source/blob/R-4-5-branch/src/main/summary.c
// R uses its native long double type when available; that is also the C++
// long double type on the supported R toolchains (including double-width
// long double on macOS arm64).
inline double reference_mean(const double* values, std::size_t size) {
  if (size == 0U) return std::numeric_limits<double>::quiet_NaN();
  const long double count = static_cast<long double>(size);
  long double sum = 0.0;
  for (std::size_t i = 0; i < size; ++i) sum += values[i];
  const bool finite_sum = std::isfinite(static_cast<double>(sum));
  if (finite_sum) {
    sum /= count;
  } else {
    sum = 0.0;
    for (std::size_t i = 0; i < size; ++i) sum += values[i] / count;
  }
  if (finite_sum && std::isfinite(static_cast<double>(sum))) {
    long double correction = 0.0;
    for (std::size_t i = 0; i < size; ++i) correction += values[i] - sum;
    sum += correction / count;
  } else if (std::isfinite(static_cast<double>(sum))) {
    long double correction = 0.0;
    for (std::size_t i = 0; i < size; ++i) correction += (values[i] - sum) / count;
    sum += correction;
  }
  return static_cast<double>(sum);
}

inline double reference_mean(const std::vector<double>& values) {
  return reference_mean(values.data(), values.size());
}

// R's prod() multiplies the double inputs in its native long-double
// accumulator and rounds once on return (rprod() in the same source file).
inline double reference_product(const std::vector<double>& values) {
  long double product = 1.0;
  for (const double value : values) product *= value;
  return static_cast<double>(product);
}

template<std::size_t Size>
inline double reference_mean(const std::array<double, Size>& values) {
  return reference_mean(values.data(), Size);
}

}  // namespace spatialess

#endif
