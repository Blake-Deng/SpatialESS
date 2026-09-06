#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

namespace {

double geometric_expression(const NumericMatrix& average,
                            const IntegerMatrix& indices,
                            int lr, int group) {
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
                         int lr, int group) {
  double factor = 1.0;
  for (int j = 0; j < indices.ncol(); ++j) {
    const int gene = indices(lr, j) - 1;
    if (gene >= 0) factor *= 1.0 + average(gene, group);
  }
  return factor;
}

double regulator_factor(const NumericMatrix& average,
                        const IntegerMatrix& indices,
                        int lr, int group, bool agonist,
                        double Kh_n, double n_power) {
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
  const double Kh_n = std::pow(Kh, n_power);
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
    receptor *= coreceptor_factor(average, co_a_indices, lr, receiver) /
      coreceptor_factor(average, co_i_indices, lr, receiver);
    if (!(ligand > 0.0) || !(receptor > 0.0)) continue;
    const double product = ligand * receptor;
    const double product_n = std::pow(product, n_power);
    double score = product_n / (Kh_n + product_n);
    if (has_agonist[lr]) {
      score *= regulator_factor(
        average, agonist_indices, lr, sender, true, Kh_n, n_power
      );
      score *= regulator_factor(
        average, agonist_indices, lr, receiver, true, Kh_n, n_power
      );
    }
    if (has_antagonist[lr]) {
      score *= regulator_factor(
        average, antagonist_indices, lr, sender, false, Kh_n, n_power
      );
      score *= regulator_factor(
        average, antagonist_indices, lr, receiver, false, Kh_n, n_power
      );
    }
    output[record] = score;
  }
  return output;
}
