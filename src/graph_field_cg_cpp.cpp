#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

double dot_product(const std::vector<double>& a, const std::vector<double>& b) {
  long double result = 0.0L;
  for (std::size_t i = 0; i < a.size(); ++i) {
    result += static_cast<long double>(a[i]) * static_cast<long double>(b[i]);
  }
  return static_cast<double>(result);
}

} // namespace

// [[Rcpp::export]]
List graph_field_cg_cpp(const IntegerVector& sender,
                        const IntegerVector& receiver,
                        const NumericVector& weight,
                        int n_cells,
                        const NumericVector& source,
                        const NumericVector& receptor,
                        double diffusion,
                        double decay,
                        double uptake,
                        double production_rate,
                        double tolerance,
                        int max_iterations) {
  const R_xlen_t edge_count = sender.size();
  if (receiver.size() != edge_count || weight.size() != edge_count) {
    stop("Graph edge vectors have inconsistent lengths.");
  }
  if (n_cells < 1 || source.size() != n_cells || receptor.size() != n_cells) {
    stop("Field vectors and graph cell count are inconsistent.");
  }

  std::vector<int> from(static_cast<std::size_t>(edge_count));
  std::vector<int> to(static_cast<std::size_t>(edge_count));
  std::vector<double> edge_weight(static_cast<std::size_t>(edge_count));
  std::vector<double> degree(static_cast<std::size_t>(n_cells), 0.0);
  for (R_xlen_t e = 0; e < edge_count; ++e) {
    const int s = sender[e] - 1;
    const int r = receiver[e] - 1;
    if (s < 0 || r < 0 || s >= n_cells || r >= n_cells) {
      stop("Graph index is outside the field vector.");
    }
    const double w = weight[e];
    if (!R_finite(w) || w < 0.0) stop("Graph weights must be finite and non-negative.");
    from[static_cast<std::size_t>(e)] = s;
    to[static_cast<std::size_t>(e)] = r;
    edge_weight[static_cast<std::size_t>(e)] = w;
    degree[static_cast<std::size_t>(s)] += w;
  }

  std::vector<double> diagonal(static_cast<std::size_t>(n_cells));
  std::vector<double> b(static_cast<std::size_t>(n_cells));
  for (int i = 0; i < n_cells; ++i) {
    if (!R_finite(source[i]) || source[i] < 0.0 ||
        !R_finite(receptor[i]) || receptor[i] < 0.0) {
      stop("Source and receptor values must be finite and non-negative.");
    }
    diagonal[static_cast<std::size_t>(i)] =
      decay + uptake * receptor[i] + diffusion * degree[static_cast<std::size_t>(i)];
    b[static_cast<std::size_t>(i)] = production_rate * source[i];
  }

  auto multiply = [&](const std::vector<double>& x, std::vector<double>& y) {
    for (int i = 0; i < n_cells; ++i) {
      y[static_cast<std::size_t>(i)] =
        diagonal[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i)];
    }
    for (R_xlen_t e = 0; e < edge_count; ++e) {
      const std::size_t ee = static_cast<std::size_t>(e);
      y[static_cast<std::size_t>(from[ee])] -=
        diffusion * edge_weight[ee] * x[static_cast<std::size_t>(to[ee])];
    }
  };

  std::vector<double> x(static_cast<std::size_t>(n_cells), 0.0);
  std::vector<double> residual = b;
  std::vector<double> z(static_cast<std::size_t>(n_cells));
  std::vector<double> direction(static_cast<std::size_t>(n_cells));
  std::vector<double> product(static_cast<std::size_t>(n_cells));
  for (int i = 0; i < n_cells; ++i) {
    z[static_cast<std::size_t>(i)] =
      residual[static_cast<std::size_t>(i)] / diagonal[static_cast<std::size_t>(i)];
    direction[static_cast<std::size_t>(i)] = z[static_cast<std::size_t>(i)];
  }

  const double b_norm = std::sqrt(dot_product(b, b));
  double residual_norm = std::sqrt(dot_product(residual, residual));
  double relative_residual = b_norm > 0.0 ? residual_norm / b_norm : 0.0;
  bool converged = relative_residual <= tolerance;
  int iterations = 0;
  double rz_old = dot_product(residual, z);

  for (int iter = 1; iter <= max_iterations && !converged; ++iter) {
    multiply(direction, product);
    const double denominator = dot_product(direction, product);
    if (!R_finite(denominator) || denominator <= 0.0) {
      stop("Conjugate-gradient breakdown: the field operator is not positive definite.");
    }
    const double step = rz_old / denominator;
    for (int i = 0; i < n_cells; ++i) {
      x[static_cast<std::size_t>(i)] += step * direction[static_cast<std::size_t>(i)];
      residual[static_cast<std::size_t>(i)] -= step * product[static_cast<std::size_t>(i)];
    }
    iterations = iter;
    residual_norm = std::sqrt(dot_product(residual, residual));
    relative_residual = b_norm > 0.0 ? residual_norm / b_norm : 0.0;
    if (relative_residual <= tolerance) {
      converged = true;
      break;
    }
    for (int i = 0; i < n_cells; ++i) {
      z[static_cast<std::size_t>(i)] =
        residual[static_cast<std::size_t>(i)] / diagonal[static_cast<std::size_t>(i)];
    }
    const double rz_new = dot_product(residual, z);
    const double beta = rz_new / rz_old;
    for (int i = 0; i < n_cells; ++i) {
      direction[static_cast<std::size_t>(i)] =
        z[static_cast<std::size_t>(i)] + beta * direction[static_cast<std::size_t>(i)];
    }
    rz_old = rz_new;
  }

  long double production_total = 0.0L;
  long double sink_total = 0.0L;
  double minimum = x.empty() ? 0.0 : x[0];
  for (int i = 0; i < n_cells; ++i) {
    production_total += static_cast<long double>(production_rate) * source[i];
    sink_total += static_cast<long double>(decay) * x[static_cast<std::size_t>(i)] +
      static_cast<long double>(uptake) * receptor[i] * x[static_cast<std::size_t>(i)];
    minimum = std::min(minimum, x[static_cast<std::size_t>(i)]);
  }
  const double production_double = static_cast<double>(production_total);
  const double sink_double = static_cast<double>(sink_total);
  const double mass_balance_relative_error =
    std::abs(production_double - sink_double) / std::max(1.0, std::abs(production_double));

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
