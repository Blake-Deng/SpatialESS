#include <Rcpp.h>
#include <chrono>
#include <vector>

#include "sparse_type7.h"

using namespace Rcpp;

// [[Rcpp::export]]
List benchmark_sparse_type7_cpp(const NumericVector& positive_values,
                                int total_count,
                                int iterations) {
  if (total_count < 1 || iterations < 1 ||
      positive_values.size() > total_count) {
    stop("Invalid sparse type-7 benchmark dimensions.");
  }
  std::vector<double> base(positive_values.begin(), positive_values.end());
  for (double value : base) {
    if (!R_finite(value) || value <= 0.0) {
      stop("positive_values must be finite and strictly positive.");
    }
  }

  using Clock = std::chrono::steady_clock;
  double sort_value = 0.0;
  double select_value = 0.0;
  double sort_seconds = 0.0;
  double select_seconds = 0.0;
  for (int iteration = 0; iteration < iterations; ++iteration) {
    std::vector<double> work = base;
    const Clock::time_point start = Clock::now();
    sort_value = spatialess::sparse_type7_trimean_sort(work, total_count);
    sort_seconds += std::chrono::duration<double>(Clock::now() - start).count();
  }
  for (int iteration = 0; iteration < iterations; ++iteration) {
    std::vector<double> work = base;
    const Clock::time_point start = Clock::now();
    select_value = spatialess::sparse_type7_trimean_select(work, total_count);
    select_seconds += std::chrono::duration<double>(Clock::now() - start).count();
  }
  return List::create(
    Named("sort_value") = sort_value,
    Named("select_value") = select_value,
    Named("sort_seconds") = sort_seconds,
    Named("select_seconds") = select_seconds
  );
}

