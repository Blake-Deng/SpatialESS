#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

struct GroupAggregate {
  double supported_edges = 0.0;
  double sum_hill = 0.0;
  double sum_spatial_score = 0.0;
};

R_xlen_t csr_offset(const NumericVector& offsets, R_xlen_t i) {
  const double value = offsets[i];
  if (!R_finite(value) || value < 0.0 || std::floor(value) != value ||
      value > static_cast<double>(R_XLEN_T_MAX)) {
    stop("Invalid CSR offset.");
  }
  return static_cast<R_xlen_t>(value);
}

} // namespace

// [[Rcpp::export]]
List aggregate_simple_lr_csr_cpp(const NumericVector& offsets,
                                 const IntegerVector& neighbors,
                                 const NumericVector& weights,
                                 const NumericVector& ligand_expression,
                                 const NumericVector& receptor_expression,
                                 const IntegerVector& group,
                                 int group_count,
                                 double Kh,
                                 double n_power,
                                 double min_expression,
                                 double max_group_pairs) {
  const R_xlen_t cell_count = ligand_expression.size();
  if (receptor_expression.size() != cell_count || group.size() != cell_count ||
      offsets.size() != cell_count + 1) {
    stop("Expression, group and CSR cell dimensions differ.");
  }
  if (weights.size() != neighbors.size() ||
      csr_offset(offsets, cell_count) != neighbors.size()) {
    stop("Invalid CSR edge vector lengths.");
  }
  if (group_count < 1) stop("group_count must be positive.");

  std::unordered_map<std::uint64_t, GroupAggregate> aggregates;
  aggregates.reserve(static_cast<std::size_t>(
    std::min<double>(max_group_pairs, 65536.0)));

  const double Kh_n = std::pow(Kh, n_power);
  double ligand_supported_cells = 0.0;
  double receptor_supported_cells = 0.0;
  double visited_out_edges = 0.0;
  double supported_edges = 0.0;

  for (R_xlen_t i = 0; i < cell_count; ++i) {
    if (receptor_expression[i] > min_expression) ++receptor_supported_cells;
  }

  for (R_xlen_t sender = 0; sender < cell_count; ++sender) {
    const double ligand = ligand_expression[sender];
    if (ligand <= min_expression) continue;
    ++ligand_supported_cells;
    const int sender_group = group[sender] - 1;
    if (sender_group < 0 || sender_group >= group_count) stop("Invalid sender group code.");
    const R_xlen_t begin = csr_offset(offsets, sender);
    const R_xlen_t end = csr_offset(offsets, sender + 1);
    visited_out_edges += static_cast<double>(end - begin);

    for (R_xlen_t p = begin; p < end; ++p) {
      const int receiver = neighbors[p] - 1;
      if (receiver < 0 || receiver >= cell_count) stop("CSR neighbor index out of range.");
      const double receptor = receptor_expression[receiver];
      if (receptor <= min_expression) continue;
      const int receiver_group = group[receiver] - 1;
      if (receiver_group < 0 || receiver_group >= group_count) {
        stop("Invalid receiver group code.");
      }
      const double molecular_product = ligand * receptor;
      const double product_n = std::pow(molecular_product, n_power);
      const double hill = product_n / (Kh_n + product_n);
      const std::uint64_t key =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(sender_group)) << 32) |
        static_cast<std::uint32_t>(receiver_group);
      auto inserted = aggregates.emplace(key, GroupAggregate{});
      GroupAggregate& aggregate = inserted.first->second;
      aggregate.supported_edges += 1.0;
      aggregate.sum_hill += hill;
      aggregate.sum_spatial_score += hill * weights[p];
      supported_edges += 1.0;
      if (static_cast<double>(aggregates.size()) > max_group_pairs) {
        stop("Sparse group-pair output exceeded max_group_pairs.");
      }
    }
  }

  std::vector<std::uint64_t> keys;
  keys.reserve(aggregates.size());
  for (const auto& item : aggregates) keys.push_back(item.first);
  std::sort(keys.begin(), keys.end());

  const R_xlen_t output_size = static_cast<R_xlen_t>(keys.size());
  IntegerVector sender_group_out(output_size);
  IntegerVector receiver_group_out(output_size);
  NumericVector edge_count_out(output_size);
  NumericVector sum_hill_out(output_size);
  NumericVector sum_score_out(output_size);
  NumericVector mean_score_out(output_size);

  for (R_xlen_t i = 0; i < output_size; ++i) {
    const std::uint64_t key = keys[static_cast<std::size_t>(i)];
    const int sender_group = static_cast<int>(key >> 32);
    const int receiver_group = static_cast<int>(key & 0xffffffffULL);
    const GroupAggregate& aggregate = aggregates.at(key);
    sender_group_out[i] = sender_group + 1;
    receiver_group_out[i] = receiver_group + 1;
    edge_count_out[i] = aggregate.supported_edges;
    sum_hill_out[i] = aggregate.sum_hill;
    sum_score_out[i] = aggregate.sum_spatial_score;
    mean_score_out[i] = aggregate.sum_spatial_score / aggregate.supported_edges;
  }

  return List::create(
    Named("group_pairs") = DataFrame::create(
      Named("sender_group") = sender_group_out,
      Named("receiver_group") = receiver_group_out,
      Named("supported_edges") = edge_count_out,
      Named("sum_hill") = sum_hill_out,
      Named("sum_spatial_score") = sum_score_out,
      Named("mean_spatial_score") = mean_score_out
    ),
    Named("ligand_supported_cells") = ligand_supported_cells,
    Named("receptor_supported_cells") = receptor_supported_cells,
    Named("visited_out_edges") = visited_out_edges,
    Named("supported_edges") = supported_edges
  );
}
