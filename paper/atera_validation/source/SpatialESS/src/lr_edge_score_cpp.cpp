#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List score_simple_lr_edge_cpp(const IntegerVector& sender,
                              const IntegerVector& receiver,
                              const NumericVector& distance,
                              const NumericVector& spatial_weight,
                              const NumericVector& ligand_expression,
                              const NumericVector& receptor_expression,
                              double Kh,
                              double n_power,
                              double min_expression,
                              double max_output_edges) {
  const R_xlen_t edge_count = sender.size();
  if (receiver.size() != edge_count || distance.size() != edge_count ||
      spatial_weight.size() != edge_count) {
    stop("Graph edge vectors have inconsistent lengths.");
  }
  if (ligand_expression.size() != receptor_expression.size()) {
    stop("Ligand and receptor expression lengths differ.");
  }
  const R_xlen_t cell_count = ligand_expression.size();
  const double Kh_n = std::pow(Kh, n_power);

  std::vector<int> out_sender;
  std::vector<int> out_receiver;
  std::vector<double> out_distance;
  std::vector<double> out_spatial_weight;
  std::vector<double> out_ligand;
  std::vector<double> out_receptor;
  std::vector<double> out_product;
  std::vector<double> out_hill;
  std::vector<double> out_score;

  const std::size_t reserve_n = static_cast<std::size_t>(
    std::min<double>(static_cast<double>(edge_count), 1000000.0));
  out_sender.reserve(reserve_n);
  out_receiver.reserve(reserve_n);
  out_distance.reserve(reserve_n);
  out_spatial_weight.reserve(reserve_n);
  out_ligand.reserve(reserve_n);
  out_receptor.reserve(reserve_n);
  out_product.reserve(reserve_n);
  out_hill.reserve(reserve_n);
  out_score.reserve(reserve_n);

  for (R_xlen_t e = 0; e < edge_count; ++e) {
    const int s = sender[e] - 1;
    const int r = receiver[e] - 1;
    if (s < 0 || r < 0 || s >= cell_count || r >= cell_count) {
      stop("Graph index is outside the expression matrix.");
    }
    const double ligand = ligand_expression[s];
    const double receptor = receptor_expression[r];
    if (ligand <= min_expression || receptor <= min_expression) continue;
    if (static_cast<double>(out_sender.size()) >= max_output_edges) {
      stop("Supported edge output exceeded max_output_edges; use aggregation or chunking.");
    }
    const double product = ligand * receptor;
    const double product_n = std::pow(product, n_power);
    const double hill = product_n / (Kh_n + product_n);
    out_sender.push_back(s + 1);
    out_receiver.push_back(r + 1);
    out_distance.push_back(distance[e]);
    out_spatial_weight.push_back(spatial_weight[e]);
    out_ligand.push_back(ligand);
    out_receptor.push_back(receptor);
    out_product.push_back(product);
    out_hill.push_back(hill);
    out_score.push_back(hill * spatial_weight[e]);
  }

  return List::create(
    Named("edges") = DataFrame::create(
      Named("sender") = wrap(out_sender),
      Named("receiver") = wrap(out_receiver),
      Named("distance") = wrap(out_distance),
      Named("spatial_weight") = wrap(out_spatial_weight),
      Named("ligand_expression") = wrap(out_ligand),
      Named("receptor_expression") = wrap(out_receptor),
      Named("molecular_product") = wrap(out_product),
      Named("hill") = wrap(out_hill),
      Named("spatial_score") = wrap(out_score)
    )
  );
}
