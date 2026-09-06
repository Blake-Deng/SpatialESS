#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

#include "sparse_type7.h"

using namespace Rcpp;

namespace {

double permutation_geometric_expression(const NumericMatrix& average,
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

double permutation_coreceptor_factor(const NumericMatrix& average,
                                     const IntegerMatrix& indices,
                                     int lr,
                                     int group) {
  double factor = 1.0;
  for (int j = 0; j < indices.ncol(); ++j) {
    const int gene = indices(lr, j) - 1;
    if (gene >= 0) factor *= 1.0 + average(gene, group);
  }
  return factor;
}

double permutation_regulator_factor(const NumericMatrix& average,
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

}  // namespace

// Return source_for_target, preserving every contiguous positive stratum.
// [[Rcpp::export]]
IntegerVector permute_within_strata_cpp(const IntegerVector& strata,
                                        int stratum_count) {
  const int cell_count = strata.size();
  if (stratum_count < 1) stop("stratum_count must be positive.");
  std::vector<std::vector<int> > members(stratum_count);
  for (int cell = 0; cell < cell_count; ++cell) {
    const int stratum = strata[cell] - 1;
    if (stratum < 0 || stratum >= stratum_count) stop("Invalid stratum code.");
    members[stratum].push_back(cell);
  }

  IntegerVector source_for_target(cell_count);
  for (int cell = 0; cell < cell_count; ++cell) source_for_target[cell] = cell + 1;
  for (int stratum = 0; stratum < stratum_count; ++stratum) {
    std::vector<int>& values = members[stratum];
    const std::vector<int> targets = values;
    for (int i = static_cast<int>(values.size()) - 1; i > 0; --i) {
      const int j = static_cast<int>(std::floor(R::unif_rand() * (i + 1)));
      std::swap(values[static_cast<std::size_t>(i)],
                values[static_cast<std::size_t>(j)]);
    }
    for (std::size_t i = 0; i < targets.size(); ++i) {
      source_for_target[targets[i]] = values[i] + 1;
    }
  }
  return source_for_target;
}

// Exact sparse triMean after assigning source profiles to fixed target groups.
// [[Rcpp::export]]
NumericMatrix group_tri_mean_dgc_permuted_cpp(
    const S4& cell_by_gene,
    const IntegerVector& target_group,
    const IntegerVector& source_for_target,
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
  if (target_group.size() != cell_count || source_for_target.size() != cell_count ||
      group_count < 1) {
    stop("Permutation, group and expression dimensions differ.");
  }

  std::vector<int> group_for_source(cell_count, -1);
  std::vector<unsigned char> source_seen(cell_count, 0);
  std::vector<int> group_sizes(group_count, 0);
  for (int target = 0; target < cell_count; ++target) {
    const int source = source_for_target[target] - 1;
    const int group = target_group[target] - 1;
    if (source < 0 || source >= cell_count || source_seen[source]) {
      stop("source_for_target must be a permutation.");
    }
    if (group < 0 || group >= group_count) stop("Invalid target group code.");
    source_seen[source] = 1;
    group_for_source[source] = group;
    ++group_sizes[group];
  }
  for (int group = 0; group < group_count; ++group) {
    if (group_sizes[group] == 0) stop("Unused target group level detected.");
  }

  NumericMatrix output(gene_count, group_count);
  std::vector<std::vector<double> > positive(group_count);
  for (int gene = 0; gene < gene_count; ++gene) {
    for (int group = 0; group < group_count; ++group) positive[group].clear();
    const int begin = column_offsets[gene];
    const int end = column_offsets[gene + 1];
    for (int p = begin; p < end; ++p) {
      const int source = row_indices[p];
      const double value = values[p];
      if (source < 0 || source >= cell_count || !R_finite(value) || value < 0.0) {
        stop("Expression must contain finite non-negative values.");
      }
      if (value > 0.0) positive[group_for_source[source]].push_back(value);
    }
    for (int group = 0; group < group_count; ++group) {
      output(gene, group) = spatialess::sparse_type7_trimean_select(
        positive[group], group_sizes[group]);
    }
  }
  return output;
}

// Matrix rows are supported group pairs; columns are LR records.
// [[Rcpp::export]]
NumericMatrix score_cellchat_group_matrix_cpp(
    const NumericMatrix& average,
    const IntegerVector& sender_group,
    const IntegerVector& receiver_group,
    const IntegerMatrix& ligand_indices,
    const IntegerMatrix& receptor_indices,
    const IntegerMatrix& co_a_indices,
    const IntegerMatrix& co_i_indices,
    const IntegerMatrix& agonist_indices,
    const IntegerMatrix& antagonist_indices,
    const LogicalVector& has_agonist,
    const LogicalVector& has_antagonist,
    double Kh,
    double n_power) {
  const int group_count = average.ncol();
  const int pair_count = sender_group.size();
  const int lr_count = ligand_indices.nrow();
  if (receiver_group.size() != pair_count || receptor_indices.nrow() != lr_count ||
      co_a_indices.nrow() != lr_count || co_i_indices.nrow() != lr_count ||
      agonist_indices.nrow() != lr_count || antagonist_indices.nrow() != lr_count ||
      has_agonist.size() != lr_count || has_antagonist.size() != lr_count) {
    stop("LR or group-support dimensions differ.");
  }
  if (!R_finite(Kh) || Kh <= 0.0 || !R_finite(n_power) || n_power <= 0.0) {
    stop("Kh and n_power must be positive and finite.");
  }

  NumericMatrix output(pair_count, lr_count);
  const double Kh_n = std::pow(Kh, n_power);
  std::vector<double> ligand(group_count), receptor(group_count);
  std::vector<double> agonist(group_count), antagonist(group_count);
  for (int lr = 0; lr < lr_count; ++lr) {
    for (int group = 0; group < group_count; ++group) {
      ligand[group] = permutation_geometric_expression(
        average, ligand_indices, lr, group);
      receptor[group] = permutation_geometric_expression(
        average, receptor_indices, lr, group);
      receptor[group] *= permutation_coreceptor_factor(
        average, co_a_indices, lr, group) /
        permutation_coreceptor_factor(average, co_i_indices, lr, group);
      agonist[group] = has_agonist[lr]
        ? permutation_regulator_factor(
            average, agonist_indices, lr, group, true, Kh_n, n_power)
        : 1.0;
      antagonist[group] = has_antagonist[lr]
        ? permutation_regulator_factor(
            average, antagonist_indices, lr, group, false, Kh_n, n_power)
        : 1.0;
    }
    for (int pair = 0; pair < pair_count; ++pair) {
      const int sender = sender_group[pair] - 1;
      const int receiver = receiver_group[pair] - 1;
      if (sender < 0 || sender >= group_count ||
          receiver < 0 || receiver >= group_count) {
        stop("Invalid supported group pair.");
      }
      const double product = ligand[sender] * receptor[receiver];
      if (!(product > 0.0)) continue;
      const double product_n = std::pow(product, n_power);
      output(pair, lr) = product_n / (Kh_n + product_n) *
        agonist[sender] * agonist[receiver] *
        antagonist[sender] * antagonist[receiver];
    }
  }
  return output;
}
