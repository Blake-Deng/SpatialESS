#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

struct CSRBucketKey {
  int sample;
  std::int64_t x;
  std::int64_t y;
  std::int64_t z;

  bool operator==(const CSRBucketKey& other) const noexcept {
    return sample == other.sample && x == other.x && y == other.y && z == other.z;
  }
};

struct CSRBucketHash {
  std::size_t operator()(const CSRBucketKey& key) const noexcept {
    std::size_t h = std::hash<int>{}(key.sample);
    auto mix = [&h](std::int64_t value) {
      const std::size_t x = std::hash<std::int64_t>{}(value);
      h ^= x + static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) + (h << 6) + (h >> 2);
    };
    mix(key.x);
    mix(key.y);
    mix(key.z);
    return h;
  }
};

double csr_weight(double distance, const std::string& method, double scale) {
  if (method == "binary") return 1.0;
  if (method == "gaussian") return std::exp(-(distance * distance) / (2.0 * scale * scale));
  if (method == "exponential") return std::exp(-distance / scale);
  throw std::invalid_argument("Unsupported weight function.");
}

R_xlen_t offset_value(const NumericVector& offsets, R_xlen_t i) {
  const double value = offsets[i];
  if (!R_finite(value) || value < 0.0 || std::floor(value) != value ||
      value > static_cast<double>(R_XLEN_T_MAX)) {
    stop("Invalid CSR offset.");
  }
  return static_cast<R_xlen_t>(value);
}

} // namespace

// [[Rcpp::export]]
List radius_graph_csr_cpp(const NumericMatrix& coords,
                          const IntegerVector& sample,
                          const IntegerVector& compartment,
                          double radius,
                          bool same_compartment,
                          std::string weight_method,
                          double weight_scale,
                          bool store_distance,
                          double max_edges) {
  const R_xlen_t n_x = coords.nrow();
  if (n_x > static_cast<R_xlen_t>(std::numeric_limits<int>::max())) {
    stop("The current CSR graph supports at most INT_MAX cells.");
  }
  const int n = static_cast<int>(n_x);
  const int dims = coords.ncol();
  if (dims < 2 || dims > 3) stop("coords must have two or three columns.");
  if (sample.size() != n || compartment.size() != n) stop("Metadata length mismatch.");
  if (!R_finite(radius) || radius <= 0.0) stop("radius must be positive and finite.");
  if (!R_finite(weight_scale) || weight_scale <= 0.0) stop("weight_scale must be positive and finite.");

  using BucketMap = std::unordered_map<CSRBucketKey, std::vector<int>, CSRBucketHash>;
  BucketMap buckets;
  buckets.reserve(static_cast<std::size_t>(n * 1.3) + 1U);
  std::vector<CSRBucketKey> cell_bucket(static_cast<std::size_t>(n));

  for (int i = 0; i < n; ++i) {
    const double x = coords(i, 0);
    const double y = coords(i, 1);
    const double z = dims == 3 ? coords(i, 2) : 0.0;
    if (!R_finite(x) || !R_finite(y) || !R_finite(z)) stop("coords contains non-finite values.");
    CSRBucketKey key{
      sample[i],
      static_cast<std::int64_t>(std::floor(x / radius)),
      static_cast<std::int64_t>(std::floor(y / radius)),
      dims == 3 ? static_cast<std::int64_t>(std::floor(z / radius)) : 0
    };
    cell_bucket[static_cast<std::size_t>(i)] = key;
    buckets[key].push_back(i);
  }

  const double radius_sq = radius * radius;
  const double numeric_tolerance =
    std::max(1.0, radius_sq) * 16.0 * std::numeric_limits<double>::epsilon();

  auto enumerate_pairs = [&](auto&& callback) {
    for (int i = 0; i < n; ++i) {
      if ((i & 65535) == 0) checkUserInterrupt();
      const CSRBucketKey base = cell_bucket[static_cast<std::size_t>(i)];
      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          const int z_min = dims == 3 ? -1 : 0;
          const int z_max = dims == 3 ? 1 : 0;
          for (int dz = z_min; dz <= z_max; ++dz) {
            const CSRBucketKey candidate{
              base.sample, base.x + dx, base.y + dy, base.z + dz
            };
            const auto found = buckets.find(candidate);
            if (found == buckets.end()) continue;
            for (const int j : found->second) {
              if (j <= i) continue;
              if (same_compartment && compartment[i] != compartment[j]) continue;
              double distance_sq = 0.0;
              for (int d = 0; d < dims; ++d) {
                const double delta = coords(i, d) - coords(j, d);
                distance_sq += delta * delta;
              }
              if (distance_sq <= radius_sq + numeric_tolerance) {
                callback(i, j, std::sqrt(std::max(0.0, distance_sq)));
              }
            }
          }
        }
      }
    }
  };

  std::vector<R_xlen_t> degree(static_cast<std::size_t>(n), 0);
  enumerate_pairs([&](int i, int j, double) {
    ++degree[static_cast<std::size_t>(i)];
    ++degree[static_cast<std::size_t>(j)];
  });

  std::vector<R_xlen_t> offsets(static_cast<std::size_t>(n) + 1U, 0);
  for (int i = 0; i < n; ++i) {
    const R_xlen_t previous = offsets[static_cast<std::size_t>(i)];
    const R_xlen_t increment = degree[static_cast<std::size_t>(i)];
    if (increment > R_XLEN_T_MAX - previous) stop("CSR edge count overflow.");
    offsets[static_cast<std::size_t>(i) + 1U] = previous + increment;
  }
  const R_xlen_t total_edges = offsets.back();
  if (static_cast<double>(total_edges) > max_edges) {
    stop("CSR edge count exceeded max_edges before allocation.");
  }

  IntegerVector neighbors(total_edges);
  NumericVector weights(total_edges);
  NumericVector distances = store_distance ? NumericVector(total_edges) : NumericVector(0);
  std::vector<R_xlen_t> cursor = offsets;
  std::vector<double> weighted_degree(static_cast<std::size_t>(n), 0.0);

  enumerate_pairs([&](int i, int j, double distance) {
    const double w = csr_weight(distance, weight_method, weight_scale);
    weighted_degree[static_cast<std::size_t>(i)] += w;
    weighted_degree[static_cast<std::size_t>(j)] += w;
    const R_xlen_t ij = cursor[static_cast<std::size_t>(i)]++;
    const R_xlen_t ji = cursor[static_cast<std::size_t>(j)]++;
    neighbors[ij] = j + 1;
    neighbors[ji] = i + 1;
    weights[ij] = w;
    weights[ji] = w;
    if (store_distance) {
      distances[ij] = distance;
      distances[ji] = distance;
    }
  });

  NumericVector offsets_r(static_cast<R_xlen_t>(offsets.size()));
  NumericVector degree_r(n);
  NumericVector weighted_degree_r(n);
  for (int i = 0; i < n; ++i) {
    offsets_r[i] = static_cast<double>(offsets[static_cast<std::size_t>(i)]);
    degree_r[i] = static_cast<double>(degree[static_cast<std::size_t>(i)]);
    weighted_degree_r[i] = weighted_degree[static_cast<std::size_t>(i)];
  }
  offsets_r[n] = static_cast<double>(total_edges);

  return List::create(
    Named("offsets") = offsets_r,
    Named("neighbors") = neighbors,
    Named("weights") = weights,
    Named("distances") = distances,
    Named("degree") = degree_r,
    Named("weighted_degree") = weighted_degree_r
  );
}

// [[Rcpp::export]]
List validate_spatial_csr_cpp(const NumericVector& offsets,
                              const IntegerVector& neighbors,
                              const NumericVector& weights,
                              const IntegerVector& sample,
                              const IntegerVector& compartment,
                              bool same_compartment) {
  const R_xlen_t n = sample.size();
  bool offsets_valid = offsets.size() == n + 1 && offset_value(offsets, 0) == 0;
  bool terminal_offset_valid = false;
  bool neighbors_in_range = true;
  bool finite_nonnegative_weights = weights.size() == neighbors.size();
  bool no_self_edges = true;
  bool no_cross_sample_edges = compartment.size() == n;
  bool compartment_barrier = compartment.size() == n;

  if (offsets_valid) {
    R_xlen_t previous = 0;
    for (R_xlen_t i = 0; i <= n; ++i) {
      const R_xlen_t current = offset_value(offsets, i);
      if (current < previous) offsets_valid = false;
      previous = current;
    }
    terminal_offset_valid = offset_value(offsets, n) == neighbors.size();
  }

  if (offsets_valid && terminal_offset_valid && weights.size() == neighbors.size()) {
    for (R_xlen_t i = 0; i < n; ++i) {
      const R_xlen_t begin = offset_value(offsets, i);
      const R_xlen_t end = offset_value(offsets, i + 1);
      for (R_xlen_t p = begin; p < end; ++p) {
        const int j = neighbors[p] - 1;
        if (j < 0 || j >= n) {
          neighbors_in_range = false;
          continue;
        }
        if (j == i) no_self_edges = false;
        if (sample[i] != sample[j]) no_cross_sample_edges = false;
        if (same_compartment && compartment[i] != compartment[j]) {
          compartment_barrier = false;
        }
        if (!R_finite(weights[p]) || weights[p] < 0.0) {
          finite_nonnegative_weights = false;
        }
      }
    }
  } else {
    neighbors_in_range = false;
    finite_nonnegative_weights = false;
    no_self_edges = false;
    no_cross_sample_edges = false;
    compartment_barrier = false;
  }

  return List::create(
    Named("offsets_valid") = offsets_valid,
    Named("terminal_offset_valid") = terminal_offset_valid,
    Named("neighbors_in_range") = neighbors_in_range,
    Named("finite_nonnegative_weights") = finite_nonnegative_weights,
    Named("no_self_edges") = no_self_edges,
    Named("no_cross_sample_edges") = no_cross_sample_edges,
    Named("compartment_barrier") = compartment_barrier
  );
}

// [[Rcpp::export]]
List materialize_csr_edges_cpp(const NumericVector& offsets,
                               const IntegerVector& neighbors,
                               const NumericVector& weights,
                               const NumericVector& distances) {
  const R_xlen_t n = offsets.size() - 1;
  const R_xlen_t edge_count = neighbors.size();
  if (weights.size() != edge_count) stop("CSR weights length mismatch.");
  const bool has_distance = distances.size() == edge_count;
  IntegerVector sender(edge_count);
  NumericVector distance_out(edge_count, NA_REAL);
  for (R_xlen_t i = 0; i < n; ++i) {
    const R_xlen_t begin = offset_value(offsets, i);
    const R_xlen_t end = offset_value(offsets, i + 1);
    for (R_xlen_t p = begin; p < end; ++p) {
      sender[p] = static_cast<int>(i) + 1;
      if (has_distance) distance_out[p] = distances[p];
    }
  }
  return List::create(
    Named("sender") = sender,
    Named("receiver") = clone(neighbors),
    Named("distance") = distance_out,
    Named("weight") = clone(weights)
  );
}

// [[Rcpp::export]]
NumericVector csr_weighted_degree_cpp(const NumericVector& offsets,
                                      const NumericVector& weights) {
  const R_xlen_t n = offsets.size() - 1;
  if (n < 0 || offset_value(offsets, n) != weights.size()) {
    stop("Invalid CSR offsets or weights.");
  }
  NumericVector degree(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    long double total = 0.0L;
    const R_xlen_t begin = offset_value(offsets, i);
    const R_xlen_t end = offset_value(offsets, i + 1);
    for (R_xlen_t p = begin; p < end; ++p) total += weights[p];
    degree[i] = static_cast<double>(total);
  }
  return degree;
}
