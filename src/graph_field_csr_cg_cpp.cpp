#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

double csr_field_dot(const std::vector<double>& a, const std::vector<double>& b) {
  long double result = 0.0L;
  for (std::size_t i = 0; i < a.size(); ++i) {
    result += static_cast<long double>(a[i]) * static_cast<long double>(b[i]);
  }
  return static_cast<double>(result);
}

R_xlen_t field_offset(const NumericVector& offsets, R_xlen_t i) {
  const double value = offsets[i];
  if (!R_finite(value) || value < 0.0 || std::floor(value) != value ||
      value > static_cast<double>(R_XLEN_T_MAX)) {
    stop("Invalid CSR offset.");
  }
  return static_cast<R_xlen_t>(value);
}

} // namespace

// [[Rcpp::export]]
List graph_field_csr_cg_cpp(const NumericVector& offsets,
                            const IntegerVector& neighbors,
                            const NumericVector& weights,
                            const NumericVector& degree,
                            const NumericVector& source,
                            const NumericVector& receptor,
                            double diffusion,
                            double decay,
                            double uptake,
                            double production_rate,
                            double tolerance,
                            int max_iterations) {
  const R_xlen_t n_x = source.size();
  if (n_x > static_cast<R_xlen_t>(std::numeric_limits<int>::max())) {
    stop("Field vector exceeds the current integer cell-index limit.");
  }
  const int n = static_cast<int>(n_x);
  if (receptor.size() != n || degree.size() != n || offsets.size() != n + 1) {
    stop("Field and CSR cell dimensions differ.");
  }
  if (weights.size() != neighbors.size() || field_offset(offsets, n) != neighbors.size()) {
    stop("Invalid CSR edge vector lengths.");
  }

  std::vector<double> diagonal(static_cast<std::size_t>(n));
  std::vector<double> b(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (!R_finite(source[i]) || source[i] < 0.0 ||
        !R_finite(receptor[i]) || receptor[i] < 0.0 ||
        !R_finite(degree[i]) || degree[i] < 0.0) {
      stop("Source, receptor and degree must be finite and non-negative.");
    }
    diagonal[static_cast<std::size_t>(i)] =
      decay + uptake * receptor[i] + diffusion * degree[i];
    b[static_cast<std::size_t>(i)] = production_rate * source[i];
  }

  auto multiply = [&](const std::vector<double>& x, std::vector<double>& y) {
    for (int i = 0; i < n; ++i) {
      y[static_cast<std::size_t>(i)] =
        diagonal[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i)];
      const R_xlen_t begin = field_offset(offsets, i);
      const R_xlen_t end = field_offset(offsets, i + 1);
      for (R_xlen_t p = begin; p < end; ++p) {
        const int j = neighbors[p] - 1;
        if (j < 0 || j >= n) stop("CSR neighbor index out of range.");
        y[static_cast<std::size_t>(i)] -=
          diffusion * weights[p] * x[static_cast<std::size_t>(j)];
      }
    }
  };

  std::vector<double> x(static_cast<std::size_t>(n), 0.0);
  std::vector<double> residual = b;
  std::vector<double> z(static_cast<std::size_t>(n));
  std::vector<double> direction(static_cast<std::size_t>(n));
  std::vector<double> product(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    z[static_cast<std::size_t>(i)] =
      residual[static_cast<std::size_t>(i)] / diagonal[static_cast<std::size_t>(i)];
    direction[static_cast<std::size_t>(i)] = z[static_cast<std::size_t>(i)];
  }

  const double b_norm = std::sqrt(csr_field_dot(b, b));
  double residual_norm = std::sqrt(csr_field_dot(residual, residual));
  double relative_residual = b_norm > 0.0 ? residual_norm / b_norm : 0.0;
  bool converged = relative_residual <= tolerance;
  int iterations = 0;
  double rz_old = csr_field_dot(residual, z);

  for (int iter = 1; iter <= max_iterations && !converged; ++iter) {
    multiply(direction, product);
    const double denominator = csr_field_dot(direction, product);
    if (!R_finite(denominator) || denominator <= 0.0) {
      stop("Conjugate-gradient breakdown: the CSR field operator is not positive definite.");
    }
    const double step = rz_old / denominator;
    for (int i = 0; i < n; ++i) {
      x[static_cast<std::size_t>(i)] += step * direction[static_cast<std::size_t>(i)];
      residual[static_cast<std::size_t>(i)] -= step * product[static_cast<std::size_t>(i)];
    }
    iterations = iter;
    residual_norm = std::sqrt(csr_field_dot(residual, residual));
    relative_residual = b_norm > 0.0 ? residual_norm / b_norm : 0.0;
    if (relative_residual <= tolerance) {
      converged = true;
      break;
    }
    for (int i = 0; i < n; ++i) {
      z[static_cast<std::size_t>(i)] =
        residual[static_cast<std::size_t>(i)] / diagonal[static_cast<std::size_t>(i)];
    }
    const double rz_new = csr_field_dot(residual, z);
    const double beta = rz_new / rz_old;
    for (int i = 0; i < n; ++i) {
      direction[static_cast<std::size_t>(i)] =
        z[static_cast<std::size_t>(i)] + beta * direction[static_cast<std::size_t>(i)];
    }
    rz_old = rz_new;
  }

  long double production_total = 0.0L;
  long double sink_total = 0.0L;
  double minimum = x.empty() ? 0.0 : x[0];
  for (int i = 0; i < n; ++i) {
    production_total += static_cast<long double>(production_rate) * source[i];
    sink_total += static_cast<long double>(decay) * x[static_cast<std::size_t>(i)] +
      static_cast<long double>(uptake) * receptor[i] * x[static_cast<std::size_t>(i)];
    minimum = std::min(minimum, x[static_cast<std::size_t>(i)]);
  }
  const double production_double = static_cast<double>(production_total);
  const double sink_double = static_cast<double>(sink_total);
  const double mass_balance_relative_error =
    std::abs(production_double - sink_double) /
    std::max(1.0, std::abs(production_double));

  return List::create(
    Named("concentration") = wrap(x),
    Named("converged") = converged,
    Named("iterations") = iterations,
    Named("residual_norm") = residual_norm,
    Named("relative_residual") = relative_residual,
    Named("minimum_concentration") = minimum,
    Named("production_total") = production_double,
    Named("sink_total") = sink_double,
    Named("mass_balance_relative_error") = mass_balance_relative_error
  );
}
