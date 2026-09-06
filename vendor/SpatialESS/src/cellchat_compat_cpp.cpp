#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "sparse_type7.h"

using namespace Rcpp;

namespace {

R_xlen_t checked_offset(const NumericVector& offsets, R_xlen_t i) {
  const double value = offsets[i];
  if (!R_finite(value) || value < 0.0 || std::floor(value) != value ||
      value > static_cast<double>(R_XLEN_T_MAX)) {
    stop("Invalid CSR offset.");
  }
  return static_cast<R_xlen_t>(value);
}

double geometric_expression(const NumericMatrix& average,
                            const IntegerMatrix& indices,
                            int lr,
                            int group) {
  double sum_log = 0.0;
  int count = 0;
  for (int j = 0; j < indices.ncol(); ++j) {
    const int gene = indices(lr, j) - 1;
    if (gene < 0) continue;
    const double value = average(gene, group);
    if (!R_finite(value) || value < 0.0) stop("Invalid group expression value.");
    if (value == 0.0) return 0.0;
    sum_log += std::log(value);
    ++count;
  }
  if (count == 0) return NA_REAL;
  return std::exp(sum_log / static_cast<double>(count));
}

double coreceptor_factor(const NumericMatrix& average,
                         const IntegerMatrix& indices,
                         int lr,
                         int group) {
  double factor = 1.0;
  for (int j = 0; j < indices.ncol(); ++j) {
    const int gene = indices(lr, j) - 1;
    if (gene < 0) continue;
    factor *= 1.0 + average(gene, group);
  }
  return factor;
}

double regulator_factor(const NumericMatrix& average,
                        const IntegerMatrix& indices,
                        int lr,
                        int group,
                        bool agonist,
                        double Kh_n,
                        double n_power) {
  double factor = 1.0;
  for (int j = 0; j < indices.ncol(); ++j) {
    const int gene = indices(lr, j) - 1;
    if (gene < 0) continue;
    const double value_n = std::pow(average(gene, group), n_power);
    factor *= agonist ? 1.0 + value_n / (Kh_n + value_n)
                      : Kh_n / (Kh_n + value_n);
  }
  return factor;
}

struct GroupEdgeAggregate {
  double count = 0.0;
  double weight = 0.0;
};

}  // namespace

// Compute CellChat's type-7 triMean without materializing sparse zeros.
// The input is a cells-by-genes dgCMatrix, so each selected gene is contiguous.
// [[Rcpp::export]]
NumericMatrix group_tri_mean_dgc_cpp(const S4& cell_by_gene,
                                     const IntegerVector& group,
                                     int group_count) {
  if (!cell_by_gene.inherits("dgCMatrix")) {
    stop("cell_by_gene must be a dgCMatrix.");
  }
  const IntegerVector dims = cell_by_gene.slot("Dim");
  const IntegerVector column_offsets = cell_by_gene.slot("p");
  const IntegerVector row_indices = cell_by_gene.slot("i");
  const NumericVector values = cell_by_gene.slot("x");
  const int cell_count = dims[0];
  const int gene_count = dims[1];
  if (group.size() != cell_count || group_count < 1) {
    stop("group and expression dimensions differ.");
  }
  if (column_offsets.size() != gene_count + 1 ||
      row_indices.size() != values.size()) {
    stop("Invalid dgCMatrix slots.");
  }

  std::vector<int> group_sizes(group_count, 0);
  for (int cell = 0; cell < cell_count; ++cell) {
    const int code = group[cell] - 1;
    if (code < 0 || code >= group_count) stop("Invalid group code.");
    ++group_sizes[code];
  }
  for (int k = 0; k < group_count; ++k) {
    if (group_sizes[k] == 0) stop("Unused group level detected.");
  }

  NumericMatrix output(gene_count, group_count);
  std::vector<std::vector<double> > positive(group_count);
  for (int gene = 0; gene < gene_count; ++gene) {
    for (int k = 0; k < group_count; ++k) positive[k].clear();
    const int begin = column_offsets[gene];
    const int end = column_offsets[gene + 1];
    for (int p = begin; p < end; ++p) {
      const int cell = row_indices[p];
      const double value = values[p];
      if (cell < 0 || cell >= cell_count || !R_finite(value) || value < 0.0) {
        stop("Expression must contain finite non-negative values.");
      }
      if (value > 0.0) positive[group[cell] - 1].push_back(value);
    }
    for (int k = 0; k < group_count; ++k) {
      output(gene, k) = spatialess::sparse_type7_trimean_select(
        positive[k], group_sizes[k]);
    }
  }
  return output;
}

// Collapse cell-level CSR edges into the spatially supported group-pair domain.
// [[Rcpp::export]]
DataFrame aggregate_group_support_csr_cpp(const NumericVector& offsets,
                                          const IntegerVector& neighbors,
                                          const NumericVector& weights,
                                          const IntegerVector& group,
                                          int group_count,
                                          double max_group_pairs) {
  const R_xlen_t cell_count = group.size();
  if (offsets.size() != cell_count + 1 || weights.size() != neighbors.size() ||
      checked_offset(offsets, cell_count) != neighbors.size()) {
    stop("Invalid CSR graph dimensions.");
  }
  std::unordered_map<std::uint64_t, GroupEdgeAggregate> aggregates;
  aggregates.reserve(static_cast<std::size_t>(
    std::min<double>(max_group_pairs, 65536.0)));

  for (R_xlen_t sender = 0; sender < cell_count; ++sender) {
    const int sender_group = group[sender] - 1;
    if (sender_group < 0 || sender_group >= group_count) stop("Invalid group code.");
    const R_xlen_t begin = checked_offset(offsets, sender);
    const R_xlen_t end = checked_offset(offsets, sender + 1);
    for (R_xlen_t p = begin; p < end; ++p) {
      const int receiver = neighbors[p] - 1;
      if (receiver < 0 || receiver >= cell_count) stop("Invalid neighbor index.");
      const int receiver_group = group[receiver] - 1;
      if (receiver_group < 0 || receiver_group >= group_count) stop("Invalid group code.");
      const std::uint64_t key =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(sender_group)) << 32) |
        static_cast<std::uint32_t>(receiver_group);
      GroupEdgeAggregate& item = aggregates[key];
      item.count += 1.0;
      item.weight += weights[p];
      if (static_cast<double>(aggregates.size()) > max_group_pairs) {
        stop("Group-support output exceeded max_group_pairs.");
      }
    }
  }

  std::vector<std::uint64_t> keys;
  keys.reserve(aggregates.size());
  for (const auto& item : aggregates) keys.push_back(item.first);
  std::sort(keys.begin(), keys.end());
  const R_xlen_t size = static_cast<R_xlen_t>(keys.size());
  IntegerVector sender_out(size), receiver_out(size);
  NumericVector count_out(size), sum_weight_out(size), mean_weight_out(size);
  for (R_xlen_t i = 0; i < size; ++i) {
    const std::uint64_t key = keys[static_cast<std::size_t>(i)];
    const GroupEdgeAggregate& item = aggregates.at(key);
    sender_out[i] = static_cast<int>(key >> 32) + 1;
    receiver_out[i] = static_cast<int>(key & 0xffffffffULL) + 1;
    count_out[i] = item.count;
    sum_weight_out[i] = item.weight;
    mean_weight_out[i] = item.weight / item.count;
  }
  return DataFrame::create(
    Named("sender_group") = sender_out,
    Named("receiver_group") = receiver_out,
    Named("supported_edges") = count_out,
    Named("sum_spatial_weight") = sum_weight_out,
    Named("mean_spatial_weight") = mean_weight_out
  );
}

// Evaluate CellChat molecular formulas only on spatially supported group pairs.
// [[Rcpp::export]]
List score_cellchat_group_support_cpp(
    const NumericMatrix& average,
    const IntegerVector& sender_group,
    const IntegerVector& receiver_group,
    const NumericVector& supported_edges,
    const NumericVector& sum_spatial_weight,
    const NumericVector& mean_spatial_weight,
    const IntegerMatrix& ligand_indices,
    const IntegerMatrix& receptor_indices,
    const IntegerMatrix& co_a_indices,
    const IntegerMatrix& co_i_indices,
    const IntegerMatrix& agonist_indices,
    const IntegerMatrix& antagonist_indices,
    const LogicalVector& has_agonist,
    const LogicalVector& has_antagonist,
    double Kh,
    double n_power,
    double max_output_records) {
  const R_xlen_t pair_count = sender_group.size();
  const int lr_count = ligand_indices.nrow();
  const int group_count = average.ncol();
  if (receiver_group.size() != pair_count || supported_edges.size() != pair_count ||
      sum_spatial_weight.size() != pair_count || mean_spatial_weight.size() != pair_count) {
    stop("Group-support columns have different lengths.");
  }
  if (receptor_indices.nrow() != lr_count || co_a_indices.nrow() != lr_count ||
      co_i_indices.nrow() != lr_count || agonist_indices.nrow() != lr_count ||
      antagonist_indices.nrow() != lr_count || has_agonist.size() != lr_count ||
      has_antagonist.size() != lr_count) {
    stop("LR component dimensions differ.");
  }
  if (!R_finite(Kh) || Kh <= 0.0 || !R_finite(n_power) || n_power <= 0.0) {
    stop("Kh and n_power must be positive and finite.");
  }

  std::vector<int> lr_out, sender_out, receiver_out;
  std::vector<double> edge_out, sum_weight_out, mean_weight_out;
  std::vector<double> molecular_out, weighted_mean_out, weighted_sum_out;
  std::vector<double> emitted(lr_count, 0.0);
  const std::size_t initial = static_cast<std::size_t>(
    std::min<double>(max_output_records, 65536.0));
  lr_out.reserve(initial);
  sender_out.reserve(initial);
  receiver_out.reserve(initial);
  edge_out.reserve(initial);
  sum_weight_out.reserve(initial);
  mean_weight_out.reserve(initial);
  molecular_out.reserve(initial);
  weighted_mean_out.reserve(initial);
  weighted_sum_out.reserve(initial);

  const double Kh_n = std::pow(Kh, n_power);
  std::vector<double> ligand(group_count), receptor(group_count);
  std::vector<double> agonist(group_count), antagonist(group_count);

  for (int lr = 0; lr < lr_count; ++lr) {
    for (int k = 0; k < group_count; ++k) {
      ligand[k] = geometric_expression(average, ligand_indices, lr, k);
      receptor[k] = geometric_expression(average, receptor_indices, lr, k);
      receptor[k] *= coreceptor_factor(average, co_a_indices, lr, k) /
                     coreceptor_factor(average, co_i_indices, lr, k);
      agonist[k] = has_agonist[lr]
        ? regulator_factor(average, agonist_indices, lr, k, true, Kh_n, n_power)
        : 1.0;
      antagonist[k] = has_antagonist[lr]
        ? regulator_factor(average, antagonist_indices, lr, k, false, Kh_n, n_power)
        : 1.0;
    }
    for (R_xlen_t p = 0; p < pair_count; ++p) {
      const int sender = sender_group[p] - 1;
      const int receiver = receiver_group[p] - 1;
      if (sender < 0 || sender >= group_count || receiver < 0 || receiver >= group_count) {
        stop("Invalid supported group pair.");
      }
      const double product = ligand[sender] * receptor[receiver];
      if (!(product > 0.0)) continue;
      const double product_n = std::pow(product, n_power);
      const double molecular = product_n / (Kh_n + product_n) *
        agonist[sender] * agonist[receiver] *
        antagonist[sender] * antagonist[receiver];
      if (!(molecular > 0.0)) continue;

      lr_out.push_back(lr + 1);
      sender_out.push_back(sender + 1);
      receiver_out.push_back(receiver + 1);
      edge_out.push_back(supported_edges[p]);
      sum_weight_out.push_back(sum_spatial_weight[p]);
      mean_weight_out.push_back(mean_spatial_weight[p]);
      molecular_out.push_back(molecular);
      weighted_mean_out.push_back(molecular * mean_spatial_weight[p]);
      weighted_sum_out.push_back(molecular * sum_spatial_weight[p]);
      emitted[lr] += 1.0;
      if (static_cast<double>(lr_out.size()) > max_output_records) {
        stop("Scored output exceeded max_output_records.");
      }
    }
  }

  return List::create(
    Named("group_pairs") = DataFrame::create(
      Named("lr_index") = wrap(lr_out),
      Named("sender_group") = wrap(sender_out),
      Named("receiver_group") = wrap(receiver_out),
      Named("supported_edges") = wrap(edge_out),
      Named("sum_spatial_weight") = wrap(sum_weight_out),
      Named("mean_spatial_weight") = wrap(mean_weight_out),
      Named("molecular_probability") = wrap(molecular_out),
      Named("mean_weighted_probability") = wrap(weighted_mean_out),
      Named("sum_weighted_probability") = wrap(weighted_sum_out)
    ),
    Named("emitted_group_pairs") = wrap(emitted)
  );
}
