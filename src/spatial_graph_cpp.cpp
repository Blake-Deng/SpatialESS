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

struct BucketKey {
  int sample;
  std::int64_t x;
  std::int64_t y;
  std::int64_t z;

  bool operator==(const BucketKey& other) const noexcept {
    return sample == other.sample && x == other.x && y == other.y && z == other.z;
  }
};

struct BucketHash {
  std::size_t operator()(const BucketKey& key) const noexcept {
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

struct Edge {
  int sender;
  int receiver;
  double distance;
  double weight;
};

double edge_weight(double distance, const std::string& method, double scale) {
  if (method == "binary") return 1.0;
  if (method == "gaussian") return std::exp(-(distance * distance) / (2.0 * scale * scale));
  if (method == "exponential") return std::exp(-distance / scale);
  throw std::invalid_argument("Unsupported weight function.");
}

} // namespace

// [[Rcpp::export]]
List radius_graph_cpp(const NumericMatrix& coords,
                      const IntegerVector& sample,
                      const IntegerVector& compartment,
                      double radius,
                      bool same_compartment,
                      bool include_self,
                      std::string weight_method,
                      double weight_scale) {
  const R_xlen_t n_x = coords.nrow();
  if (n_x > static_cast<R_xlen_t>(std::numeric_limits<int>::max())) {
    stop("The current in-memory graph prototype supports at most INT_MAX cells.");
  }
  const int n = static_cast<int>(n_x);
  const int dims = coords.ncol();
  if (dims < 2 || dims > 3) stop("coords must have two or three columns.");
  if (sample.size() != n || compartment.size() != n) stop("Metadata length mismatch.");
  if (!R_finite(radius) || radius <= 0.0) stop("radius must be positive and finite.");
  if (!R_finite(weight_scale) || weight_scale <= 0.0) stop("weight_scale must be positive and finite.");

  using BucketMap = std::unordered_map<BucketKey, std::vector<int>, BucketHash>;
  BucketMap buckets;
  buckets.reserve(static_cast<std::size_t>(n * 1.3) + 1U);
  std::vector<BucketKey> cell_bucket(static_cast<std::size_t>(n));

  for (int i = 0; i < n; ++i) {
    const double x = coords(i, 0);
    const double y = coords(i, 1);
    const double z = dims == 3 ? coords(i, 2) : 0.0;
    if (!R_finite(x) || !R_finite(y) || !R_finite(z)) stop("coords contains non-finite values.");
    BucketKey key{
      sample[i],
      static_cast<std::int64_t>(std::floor(x / radius)),
      static_cast<std::int64_t>(std::floor(y / radius)),
      dims == 3 ? static_cast<std::int64_t>(std::floor(z / radius)) : 0
    };
    cell_bucket[static_cast<std::size_t>(i)] = key;
    buckets[key].push_back(i);
  }

  std::vector<Edge> edges;
  const double radius_sq = radius * radius;
  const double tol = std::max(1.0, radius_sq) * 16.0 * std::numeric_limits<double>::epsilon();

  if (include_self) {
    edges.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) edges.push_back(Edge{i, i, 0.0, edge_weight(0.0, weight_method, weight_scale)});
  }

  for (int i = 0; i < n; ++i) {
    const BucketKey base = cell_bucket[static_cast<std::size_t>(i)];
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        const int z_min = dims == 3 ? -1 : 0;
        const int z_max = dims == 3 ? 1 : 0;
        for (int dz = z_min; dz <= z_max; ++dz) {
          const BucketKey candidate{base.sample, base.x + dx, base.y + dy, base.z + dz};
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
            if (distance_sq > radius_sq + tol) continue;
            const double distance = std::sqrt(std::max(0.0, distance_sq));
            const double w = edge_weight(distance, weight_method, weight_scale);
            edges.push_back(Edge{i, j, distance, w});
            edges.push_back(Edge{j, i, distance, w});
          }
        }
      }
    }
  }

  std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
    if (a.sender != b.sender) return a.sender < b.sender;
    return a.receiver < b.receiver;
  });

  const R_xlen_t m = static_cast<R_xlen_t>(edges.size());
  IntegerVector sender(m), receiver(m);
  NumericVector distance(m), weight(m);
  for (R_xlen_t e = 0; e < m; ++e) {
    const Edge& edge = edges[static_cast<std::size_t>(e)];
    sender[e] = edge.sender + 1;
    receiver[e] = edge.receiver + 1;
    distance[e] = edge.distance;
    weight[e] = edge.weight;
  }
  return List::create(
    Named("sender") = sender,
    Named("receiver") = receiver,
    Named("distance") = distance,
    Named("weight") = weight
  );
}

