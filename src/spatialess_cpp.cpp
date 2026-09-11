#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

struct DGCView {
  IntegerVector dims;
  IntegerVector p;
  IntegerVector i;
  NumericVector x;

  explicit DGCView(const S4& matrix)
    : dims(matrix.slot("Dim")), p(matrix.slot("p")),
      i(matrix.slot("i")), x(matrix.slot("x")) {
    if (!matrix.inherits("dgCMatrix")) stop("cell_by_gene must be a dgCMatrix.");
    if (dims.size() != 2 || p.size() != dims[1] + 1 || i.size() != x.size()) {
      stop("Invalid dgCMatrix slots.");
    }
  }

  int cells() const { return dims[0]; }
  int genes() const { return dims[1]; }
};

R_xlen_t csr_offset(const NumericVector& offsets, R_xlen_t index) {
  const double value = offsets[index];
  if (!R_finite(value) || value < 0.0 || std::floor(value) != value ||
      value > static_cast<double>(R_XLEN_T_MAX)) {
    stop("Invalid CSR offset.");
  }
  return static_cast<R_xlen_t>(value);
}

std::uint64_t group_key(int sender, int receiver) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(sender)) << 32) |
    static_cast<std::uint32_t>(receiver);
}

int key_sender(std::uint64_t key) {
  return static_cast<int>(key >> 32);
}

int key_receiver(std::uint64_t key) {
  return static_cast<int>(key & 0xffffffffULL);
}

std::vector<int> component_genes(const IntegerMatrix& indices, int row,
                                 int gene_count) {
  std::vector<int> genes;
  genes.reserve(indices.ncol());
  for (int column = 0; column < indices.ncol(); ++column) {
    const int gene = indices(row, column) - 1;
    if (gene < 0) continue;
    if (gene >= gene_count) stop("LR component gene index is out of range.");
    genes.push_back(gene);
  }
  return genes;
}

void fill_geometric_expression(const DGCView& matrix,
                               const IntegerMatrix& indices,
                               int row,
                               std::vector<double>& output,
                               std::vector<double>& log_sum,
                               std::vector<int>& present) {
  const std::vector<int> genes = component_genes(indices, row, matrix.genes());
  if (genes.empty()) stop("Ligand/receptor component has no genes.");
  std::fill(output.begin(), output.end(), 0.0);
  if (genes.size() == 1U) {
    const int gene = genes.front();
    for (int q = matrix.p[gene]; q < matrix.p[gene + 1]; ++q) {
      const double value = matrix.x[q];
      if (!R_finite(value) || value < 0.0) stop("Invalid expression value.");
      if (value > 0.0) output[static_cast<std::size_t>(matrix.i[q])] = value;
    }
    return;
  }

  std::fill(log_sum.begin(), log_sum.end(), 0.0);
  std::fill(present.begin(), present.end(), 0);
  for (const int gene : genes) {
    for (int q = matrix.p[gene]; q < matrix.p[gene + 1]; ++q) {
      const int cell = matrix.i[q];
      const double value = matrix.x[q];
      if (!R_finite(value) || value < 0.0) stop("Invalid expression value.");
      if (value > 0.0) {
        log_sum[static_cast<std::size_t>(cell)] += std::log(value);
        ++present[static_cast<std::size_t>(cell)];
      }
    }
  }
  const int required = static_cast<int>(genes.size());
  for (int cell = 0; cell < matrix.cells(); ++cell) {
    if (present[static_cast<std::size_t>(cell)] == required) {
      output[static_cast<std::size_t>(cell)] = std::exp(
        log_sum[static_cast<std::size_t>(cell)] / static_cast<double>(required));
    }
  }
}

enum class FactorMode { Coreceptor, Agonist, Antagonist };

void fill_factor(const DGCView& matrix,
                 const IntegerMatrix& indices,
                 int row,
                 FactorMode mode,
                 double Kh_n,
                 double n_power,
                 std::vector<double>& output) {
  std::fill(output.begin(), output.end(), 1.0);
  const std::vector<int> genes = component_genes(indices, row, matrix.genes());
  for (const int gene : genes) {
    for (int q = matrix.p[gene]; q < matrix.p[gene + 1]; ++q) {
      const int cell = matrix.i[q];
      const double value = matrix.x[q];
      if (!R_finite(value) || value < 0.0) stop("Invalid expression value.");
      double factor = 1.0;
      if (mode == FactorMode::Coreceptor) {
        factor = 1.0 + value;
      } else {
        const double value_n = std::pow(value, n_power);
        factor = mode == FactorMode::Agonist
          ? 1.0 + value_n / (Kh_n + value_n)
          : Kh_n / (Kh_n + value_n);
      }
      output[static_cast<std::size_t>(cell)] *= factor;
    }
  }
}

struct Aggregate {
  double sum = 0.0;
  double count = 0.0;
};

struct ActiveEdge {
  int sender;
  int receiver;
  double probability;
};

void add_edge(std::unordered_map<std::uint64_t, Aggregate>& aggregate,
              int sender_group,
              int receiver_group,
              double probability) {
  Aggregate& item = aggregate[group_key(sender_group, receiver_group)];
  item.sum += probability;
  item.count += 1.0;
}

LogicalMatrix v3_percent_mask(const std::vector<int>& ligand_positive,
                             const std::vector<int>& receptor_positive,
                             const std::vector<int>& group_sizes,
                             double min_percent,
                             const Function& filter) {
  const int groups = static_cast<int>(group_sizes.size());
  LogicalMatrix mask(groups, 2);
  // All formatted non-negative fractions pass the zero threshold in V3.
  if (min_percent == 0.0) {
    std::fill(mask.begin(), mask.end(), true);
    return mask;
  }
  NumericMatrix positive(groups, 2);
  for (int g = 0; g < groups; ++g) {
    positive(g, 0) = ligand_positive[static_cast<std::size_t>(g)];
    positive(g, 1) = receptor_positive[static_cast<std::size_t>(g)];
  }
  // Format only the small group-level table, never a cell-by-cell matrix.
  return filter(positive, wrap(group_sizes), min_percent);
}

}  // namespace

// [[Rcpp::export]]
List spatialess_stream_cpp(
    const S4& cell_by_gene,
    const NumericVector& offsets,
    const IntegerVector& neighbors,
    const NumericVector& distances,
    const IntegerVector& group,
    int group_count,
    const IntegerMatrix& ligand_indices,
    const IntegerMatrix& receptor_indices,
    const IntegerMatrix& co_a_indices,
    const IntegerMatrix& co_i_indices,
    const IntegerMatrix& agonist_indices,
    const IntegerMatrix& antagonist_indices,
    const LogicalVector& has_agonist,
    const LogicalVector& has_antagonist,
    const LogicalVector& contact_lr,
    const IntegerMatrix& permutations,
    double Kh,
    double n_power,
    double scale_distance,
    double contact_threshold,
    double min_percent,
    int min_cells_sr,
    bool use_agan,
    double max_output_records,
    double max_active_edges_per_lr) {
  const DGCView expression(cell_by_gene);
  const int cell_count = expression.cells();
  const int lr_count = ligand_indices.nrow();
  if (group.size() != cell_count || offsets.size() != cell_count + 1 ||
      csr_offset(offsets, cell_count) != neighbors.size() ||
      distances.size() != neighbors.size()) {
    stop("Expression, group and CSR graph dimensions differ.");
  }
  if (group_count < 1 || receptor_indices.nrow() != lr_count ||
      co_a_indices.nrow() != lr_count || co_i_indices.nrow() != lr_count ||
      agonist_indices.nrow() != lr_count || antagonist_indices.nrow() != lr_count ||
      has_agonist.size() != lr_count || has_antagonist.size() != lr_count ||
      contact_lr.size() != lr_count) {
    stop("Invalid LR component dimensions.");
  }
  if (permutations.nrow() != cell_count) {
    stop("permutations must have one row per cell.");
  }
  if (!R_finite(Kh) || Kh <= 0.0 || !R_finite(n_power) || n_power <= 0.0 ||
      !R_finite(scale_distance) || scale_distance <= 0.0 ||
      !R_finite(contact_threshold) || contact_threshold < 0.0 ||
      !R_finite(min_percent) || min_percent < 0.0 || min_percent > 1.0 ||
      min_cells_sr < 1) {
    stop("Invalid probability or filtering parameter.");
  }
  for (int cell = 0; cell < cell_count; ++cell) {
    if (group[cell] < 1 || group[cell] > group_count) stop("Invalid group code.");
  }

  double minimum_distance = std::numeric_limits<double>::infinity();
  for (R_xlen_t edge = 0; edge < distances.size(); ++edge) {
    const double distance = distances[edge];
    if (!R_finite(distance) || distance <= 0.0) stop("CSR distances must be positive.");
    minimum_distance = std::min(minimum_distance, distance);
  }
  if (!R_finite(minimum_distance)) stop("The spatial graph has no non-self edges.");
  const double self_spatial_weight = 1.0 / (minimum_distance * scale_distance);
  const double Kh_n = std::pow(Kh, n_power);
  const int nboot = permutations.ncol();
  const Environment package = Environment::namespace_env("SpatialESS");
  const Function percent_filter = package[".spatialess_v3_percent_mask"];

  std::vector<int> lr_out, sender_group_out, receiver_group_out;
  std::vector<double> probability_out, pvalue_out, edge_count_out;
  std::vector<double> active_edges_per_lr(lr_count, 0.0);
  std::vector<double> ligand(cell_count), receptor(cell_count);
  std::vector<double> co_a(cell_count), co_i(cell_count);
  std::vector<double> agonist(cell_count), antagonist(cell_count);
  std::vector<double> log_sum(cell_count);
  std::vector<int> present(cell_count);
  std::vector<double> sender_links(cell_count), receiver_links(cell_count);
  std::vector<int> group_sizes(group_count), ligand_positive(group_count),
    receptor_positive(group_count);
  std::vector<double> sender_group_links(group_count), receiver_group_links(group_count);
  std::vector<int> boot_group(cell_count), boot_group_sizes(group_count),
    boot_ligand_positive(group_count), boot_receptor_positive(group_count);
  std::vector<double> boot_sender_links(group_count), boot_receiver_links(group_count);

  for (int lr = 0; lr < lr_count; ++lr) {
    if ((lr & 15) == 0) checkUserInterrupt();
    fill_geometric_expression(expression, ligand_indices, lr, ligand, log_sum, present);
    fill_geometric_expression(expression, receptor_indices, lr, receptor, log_sum, present);
    fill_factor(expression, co_a_indices, lr, FactorMode::Coreceptor,
                Kh_n, n_power, co_a);
    fill_factor(expression, co_i_indices, lr, FactorMode::Coreceptor,
                Kh_n, n_power, co_i);
    for (int cell = 0; cell < cell_count; ++cell) {
      receptor[static_cast<std::size_t>(cell)] *=
        co_a[static_cast<std::size_t>(cell)] / co_i[static_cast<std::size_t>(cell)];
    }
    if (use_agan && has_agonist[lr]) {
      fill_factor(expression, agonist_indices, lr, FactorMode::Agonist,
                  Kh_n, n_power, agonist);
    } else {
      std::fill(agonist.begin(), agonist.end(), 1.0);
    }
    if (use_agan && has_antagonist[lr]) {
      fill_factor(expression, antagonist_indices, lr, FactorMode::Antagonist,
                  Kh_n, n_power, antagonist);
    } else {
      std::fill(antagonist.begin(), antagonist.end(), 1.0);
    }

    std::fill(sender_links.begin(), sender_links.end(), 0.0);
    std::fill(receiver_links.begin(), receiver_links.end(), 0.0);
    std::fill(group_sizes.begin(), group_sizes.end(), 0);
    std::fill(ligand_positive.begin(), ligand_positive.end(), 0);
    std::fill(receptor_positive.begin(), receptor_positive.end(), 0);
    std::unordered_map<std::uint64_t, Aggregate> observed;
    observed.reserve(1024);
    std::vector<ActiveEdge> active_edges;
    if (nboot > 0) active_edges.reserve(
      static_cast<std::size_t>(std::min<double>(max_active_edges_per_lr, 65536.0)));

    for (int cell = 0; cell < cell_count; ++cell) {
      const int code = group[cell] - 1;
      ++group_sizes[static_cast<std::size_t>(code)];
      if (ligand[static_cast<std::size_t>(cell)] > 0.0) {
        ++ligand_positive[static_cast<std::size_t>(code)];
      }
      if (receptor[static_cast<std::size_t>(cell)] > 0.0) {
        ++receptor_positive[static_cast<std::size_t>(code)];
      }
    }

    auto process_edge = [&](int sender, int receiver, double spatial_weight) {
      const double product = ligand[static_cast<std::size_t>(sender)] *
        receptor[static_cast<std::size_t>(receiver)];
      if (!(product > 0.0)) return;
      const double product_n = std::pow(product, n_power);
      double probability = product_n / (Kh_n + product_n) * spatial_weight;
      if (use_agan) {
        probability *= agonist[static_cast<std::size_t>(sender)] *
          agonist[static_cast<std::size_t>(receiver)] *
          antagonist[static_cast<std::size_t>(sender)] *
          antagonist[static_cast<std::size_t>(receiver)];
      }
      if (!(probability > 0.0) || !R_finite(probability)) return;
      add_edge(observed, group[sender] - 1, group[receiver] - 1, probability);
      sender_links[static_cast<std::size_t>(sender)] += 1.0;
      receiver_links[static_cast<std::size_t>(receiver)] += 1.0;
      if (nboot > 0) {
        active_edges.push_back(ActiveEdge{sender, receiver, probability});
        if (static_cast<double>(active_edges.size()) > max_active_edges_per_lr) {
          stop("Active edge count exceeded max_active_edges_per_lr.");
        }
      }
    };

    for (int sender = 0; sender < cell_count; ++sender) {
      if (!(ligand[static_cast<std::size_t>(sender)] > 0.0)) continue;
      const R_xlen_t begin = csr_offset(offsets, sender);
      const R_xlen_t end = csr_offset(offsets, sender + 1);
      for (R_xlen_t edge = begin; edge < end; ++edge) {
        const int receiver = neighbors[edge] - 1;
        if (receiver < 0 || receiver >= cell_count) stop("Invalid neighbor index.");
        const double distance = distances[edge];
        if (contact_lr[lr] && distance > contact_threshold) continue;
        process_edge(sender, receiver, 1.0 / (distance * scale_distance));
      }
      process_edge(sender, sender, self_spatial_weight);
    }
    if (nboot > 0) {
      std::sort(active_edges.begin(), active_edges.end(),
        [](const ActiveEdge& left, const ActiveEdge& right) {
          if (left.sender != right.sender) return left.sender < right.sender;
          return left.receiver < right.receiver;
        });
      observed.clear();
      observed.reserve(1024);
      for (const ActiveEdge& edge : active_edges) {
        add_edge(observed, group[edge.sender] - 1,
                 group[edge.receiver] - 1, edge.probability);
      }
    }
    active_edges_per_lr[static_cast<std::size_t>(lr)] =
      nboot > 0 ? static_cast<double>(active_edges.size()) :
      std::accumulate(observed.begin(), observed.end(), 0.0,
        [](double total, const auto& item) { return total + item.second.count; });

    std::fill(sender_group_links.begin(), sender_group_links.end(), 0.0);
    std::fill(receiver_group_links.begin(), receiver_group_links.end(), 0.0);
    for (int cell = 0; cell < cell_count; ++cell) {
      const int code = group[cell] - 1;
      sender_group_links[static_cast<std::size_t>(code)] +=
        sender_links[static_cast<std::size_t>(cell)];
      receiver_group_links[static_cast<std::size_t>(code)] +=
        receiver_links[static_cast<std::size_t>(cell)];
    }

    std::vector<std::uint64_t> keys;
    keys.reserve(observed.size());
    LogicalMatrix observed_percent(group_count, 2);
    if (!observed.empty()) {
      observed_percent = v3_percent_mask(ligand_positive, receptor_positive,
                                        group_sizes, min_percent, percent_filter);
    }
    for (const auto& item : observed) {
      const int sender_group = key_sender(item.first);
      const int receiver_group = key_receiver(item.first);
      const bool percent_ok = observed_percent(sender_group, 0) &&
        observed_percent(receiver_group, 1);
      const bool cells_ok =
        sender_group_links[static_cast<std::size_t>(sender_group)] >= min_cells_sr &&
        receiver_group_links[static_cast<std::size_t>(receiver_group)] >= min_cells_sr;
      if (percent_ok && cells_ok && item.second.count > 0.0) keys.push_back(item.first);
    }
    std::sort(keys.begin(), keys.end());
    std::vector<int> reject(keys.size(), 0);

    for (int boot = 0; boot < nboot; ++boot) {
      std::fill(boot_group_sizes.begin(), boot_group_sizes.end(), 0);
      std::fill(boot_ligand_positive.begin(), boot_ligand_positive.end(), 0);
      std::fill(boot_receptor_positive.begin(), boot_receptor_positive.end(), 0);
      std::fill(boot_sender_links.begin(), boot_sender_links.end(), 0.0);
      std::fill(boot_receiver_links.begin(), boot_receiver_links.end(), 0.0);
      for (int cell = 0; cell < cell_count; ++cell) {
        const int source = permutations(cell, boot) - 1;
        if (source < 0 || source >= cell_count) stop("Invalid permutation index.");
        const int code = group[source] - 1;
        boot_group[static_cast<std::size_t>(cell)] = code;
        ++boot_group_sizes[static_cast<std::size_t>(code)];
        if (ligand[static_cast<std::size_t>(cell)] > 0.0) {
          ++boot_ligand_positive[static_cast<std::size_t>(code)];
        }
        if (receptor[static_cast<std::size_t>(cell)] > 0.0) {
          ++boot_receptor_positive[static_cast<std::size_t>(code)];
        }
        boot_sender_links[static_cast<std::size_t>(code)] +=
          sender_links[static_cast<std::size_t>(cell)];
        boot_receiver_links[static_cast<std::size_t>(code)] +=
          receiver_links[static_cast<std::size_t>(cell)];
      }
      std::unordered_map<std::uint64_t, Aggregate> boot_aggregate;
      boot_aggregate.reserve(observed.size());
      for (const ActiveEdge& edge : active_edges) {
        add_edge(boot_aggregate,
                 boot_group[static_cast<std::size_t>(edge.sender)],
                 boot_group[static_cast<std::size_t>(edge.receiver)],
                 edge.probability);
      }
      LogicalMatrix boot_percent(group_count, 2);
      if (!keys.empty()) {
        boot_percent = v3_percent_mask(boot_ligand_positive, boot_receptor_positive,
                                      boot_group_sizes, min_percent, percent_filter);
      }
      for (std::size_t index = 0; index < keys.size(); ++index) {
        const std::uint64_t key = keys[index];
        const int sender_group = key_sender(key);
        const int receiver_group = key_receiver(key);
        const bool percent_ok = boot_percent(sender_group, 0) &&
          boot_percent(receiver_group, 1);
        const bool cells_ok =
          boot_sender_links[static_cast<std::size_t>(sender_group)] >= min_cells_sr &&
          boot_receiver_links[static_cast<std::size_t>(receiver_group)] >= min_cells_sr;
        double boot_probability = 0.0;
        const auto found = boot_aggregate.find(key);
        if (percent_ok && cells_ok && found != boot_aggregate.end() &&
            found->second.count > 0.0) {
          boot_probability = found->second.sum / found->second.count;
        }
        const Aggregate& observed_item = observed.at(key);
        const double observed_probability = observed_item.sum / observed_item.count;
        if (boot_probability > observed_probability) ++reject[index];
      }
    }

    for (std::size_t index = 0; index < keys.size(); ++index) {
      const std::uint64_t key = keys[index];
      const Aggregate& item = observed.at(key);
      lr_out.push_back(lr + 1);
      sender_group_out.push_back(key_sender(key) + 1);
      receiver_group_out.push_back(key_receiver(key) + 1);
      probability_out.push_back(item.sum / item.count);
      edge_count_out.push_back(item.count);
      pvalue_out.push_back(nboot > 0
        ? static_cast<double>(reject[index]) / static_cast<double>(nboot)
        : NA_REAL);
      if (static_cast<double>(lr_out.size()) > max_output_records) {
        stop("Output exceeded max_output_records.");
      }
    }
  }

  return List::create(
    Named("records") = DataFrame::create(
      Named("lr_index") = wrap(lr_out),
      Named("sender_group") = wrap(sender_group_out),
      Named("receiver_group") = wrap(receiver_group_out),
      Named("probability") = wrap(probability_out),
      Named("pvalue") = wrap(pvalue_out),
      Named("active_edges") = wrap(edge_count_out)
    ),
    Named("active_edges_per_lr") = wrap(active_edges_per_lr),
    Named("minimum_nonzero_distance") = minimum_distance,
    Named("self_spatial_weight") = self_spatial_weight
  );
}
