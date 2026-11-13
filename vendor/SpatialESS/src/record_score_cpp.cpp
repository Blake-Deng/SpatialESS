#include <Rcpp.h>
#include <cmath>

#include "reference_arithmetic.h"

using namespace Rcpp;

namespace {

double geometric_expression(const NumericMatrix& average,
                            const IntegerMatrix& indices,
                            int lr, int group) {
  std::vector<double> log_values;
  log_values.reserve(indices.ncol());
  double single_value = 0.0;
  int count = 0;
  for (int j = 0; j < indices.ncol(); ++j) {
    const int gene = indices(lr, j) - 1;
    if (gene < 0) continue;
    const double value = average(gene, group);
    if (!R_finite(value) || value < 0.0) stop("Invalid group expression value.");
    if (value == 0.0) return 0.0;
    single_value = value;
    log_values.push_back(std::log(value));
    ++count;
  }
  if (count == 0) return NA_REAL;
  if (count == 1) {
    bool direct_gene = true;
    if (indices.hasAttribute("spatialess_single_gene")) {
      const LogicalVector single_gene = indices.attr("spatialess_single_gene");
      if (single_gene.size() != indices.nrow() || single_gene[lr] == NA_LOGICAL) {
        stop("Invalid single-gene component metadata.");
      }
      direct_gene = single_gene[lr];
    }
    if (direct_gene) return single_value;
  }
  return std::exp(spatialess::reference_mean(log_values));
}

double coreceptor_factor(const NumericMatrix& average,
                         const IntegerMatrix& indices,
                         int lr, int group) {
  std::vector<double> factors;
  factors.reserve(indices.ncol());
  for (int j = 0; j < indices.ncol(); ++j) {
    const int gene = indices(lr, j) - 1;
    if (gene >= 0) factors.push_back(1.0 + average(gene, group));
  }
  return spatialess::reference_product(factors);
}

double regulator_factor(const NumericMatrix& average,
                        const IntegerMatrix& indices,
                        int lr, int group, bool agonist,
                        double Kh_n, double n_power) {
  std::vector<double> factors;
  factors.reserve(indices.ncol());
  for (int j = 0; j < indices.ncol(); ++j) {
    const int gene = indices(lr, j) - 1;
    if (gene < 0) continue;
    const double value_n = R_pow(average(gene, group), n_power);
    factors.push_back(agonist ? 1.0 + value_n / (Kh_n + value_n)
                             : Kh_n / (Kh_n + value_n));
  }
  return spatialess::reference_product(factors);
}

}  // namespace

// Score only the requested support-pair/LR records.
// [[Rcpp::export]]
NumericVector score_cellchat_group_records_cpp(
    const NumericMatrix& average,
    const IntegerVector& sender_group,
    const IntegerVector& receiver_group,
    const IntegerVector& record_lr,
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
  const int record_count = sender_group.size();
  const int group_count = average.ncol();
  const int lr_count = ligand_indices.nrow();
  if (receiver_group.size() != record_count || record_lr.size() != record_count ||
      receptor_indices.nrow() != lr_count || co_a_indices.nrow() != lr_count ||
      co_i_indices.nrow() != lr_count || agonist_indices.nrow() != lr_count ||
      antagonist_indices.nrow() != lr_count || has_agonist.size() != lr_count ||
      has_antagonist.size() != lr_count) {
    stop("Record, LR or component dimensions differ.");
  }
  if (!R_finite(Kh) || Kh <= 0.0 || !R_finite(n_power) || n_power <= 0.0) {
    stop("Kh and n_power must be positive and finite.");
  }
  NumericVector output(record_count);
  const double Kh_n = R_pow(Kh, n_power);
  for (int record = 0; record < record_count; ++record) {
    const int sender = sender_group[record] - 1;
    const int receiver = receiver_group[record] - 1;
    const int lr = record_lr[record] - 1;
    if (sender < 0 || sender >= group_count || receiver < 0 ||
        receiver >= group_count || lr < 0 || lr >= lr_count) {
      stop("Invalid group or LR record index.");
    }
    const double ligand = geometric_expression(
      average, ligand_indices, lr, sender
    );
    double receptor = geometric_expression(
      average, receptor_indices, lr, receiver
    );
    receptor = (receptor * coreceptor_factor(
      average, co_a_indices, lr, receiver)) /
      coreceptor_factor(average, co_i_indices, lr, receiver);
    if (!(ligand > 0.0) || !(receptor > 0.0)) continue;
    const double product = ligand * receptor;
    const double product_n = R_pow(product, n_power);
    double score = product_n / (Kh_n + product_n);
    if (has_agonist[lr]) {
      const double sender_factor = regulator_factor(
        average, agonist_indices, lr, sender, true, Kh_n, n_power);
      const double receiver_factor = regulator_factor(
        average, agonist_indices, lr, receiver, true, Kh_n, n_power);
      score *= sender_factor * receiver_factor;
    }
    if (has_antagonist[lr]) {
      const double sender_factor = regulator_factor(
        average, antagonist_indices, lr, sender, false, Kh_n, n_power);
      const double receiver_factor = regulator_factor(
        average, antagonist_indices, lr, receiver, false, Kh_n, n_power);
      score *= sender_factor * receiver_factor;
    }
    output[record] = score;
  }
  return output;
}
