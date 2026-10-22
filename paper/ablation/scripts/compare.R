args <- commandArgs(TRUE)
root <- normalizePath(args[1]); phase <- args[2]; id <- args[3]; nboot <- as.integer(args[4])
target_variant <- args[5]; repetition <- as.integer(args[6])
reference_variant <- if (length(args) > 6L) args[7] else 'full'
key <- function(x) paste(x$sender_group_name, x$receiver_group_name,
                        x$interaction_name, sep = '\037')
directory <- function(variant, rep) file.path(root, 'results', phase,
  sprintf('%s_%s_b%d_rep%d', id, variant, nboot, rep))
ref_dir <- directory(reference_variant, if (phase %in% c('tiny', 'preflight')) 0L else repetition)
target_dir <- directory(target_variant, repetition)
ref <- readRDS(file.path(ref_dir, 'records.rds'))
target <- readRDS(file.path(target_dir, 'records.rds'))
rk <- key(ref); tk <- key(target)
stopifnot(!anyDuplicated(rk), !anyDuplicated(tk))
keys_equal <- setequal(rk, tk)
common <- intersect(rk, tk)
r <- ref[match(common, rk), ]; t <- target[match(common, tk), ]
prob_error <- if (length(common)) max(abs(r$probability - t$probability)) else 0
p_error <- if (length(common)) max(abs(r$pvalue - t$pvalue)) else 0
p_equal <- identical(r$pvalue, t$pvalue)
disagreement <- sum((r$pvalue < .05) != (t$pvalue < .05))
disagreement_le <- sum((r$pvalue <= .05) != (t$pvalue <= .05))
edges_equal <- identical(r$active_edges, t$active_edges)
rs <- read.delim(file.path(ref_dir, 'summary.tsv'))
ts <- read.delim(file.path(target_dir, 'summary.tsv'))
graph_equal <- rs$directed_edges == ts$directed_edges
graph_max_error <- NA_real_
if (file.exists(file.path(ref_dir, 'graph.rds')) && file.exists(file.path(target_dir, 'graph.rds'))) {
  edges <- function(g) {
    if (!is.null(g$dense_distances)) {
      index <- which(g$dense_distances > 0, arr.ind = TRUE)
      d <- data.frame(sender = index[, 1], receiver = index[, 2],
                      distance = g$dense_distances[index])
    } else {
      sender <- if (is.null(g$senders)) rep(seq_len(length(g$offsets) - 1L), diff(g$offsets)) else g$senders
      d <- data.frame(sender = sender, receiver = g$neighbors, distance = g$distances)
    }
    d[order(d$sender, d$receiver), ]
  }
  re <- edges(readRDS(file.path(ref_dir, 'graph.rds')))
  te <- edges(readRDS(file.path(target_dir, 'graph.rds')))
  graph_equal <- graph_equal && identical(re$sender, te$sender) && identical(re$receiver, te$receiver)
  graph_max_error <- if (graph_equal && nrow(re)) max(abs(re$distance - te$distance)) else if (graph_equal) 0 else Inf
  graph_equal <- graph_equal && graph_max_error == 0
}
pass <- keys_equal && graph_equal && edges_equal && prob_error <= 1e-15 && p_equal && disagreement == 0 && disagreement_le == 0
row <- data.frame(case = id, phase = phase, variant = target_variant, repetition = repetition,
  nboot = nboot, reference_variant = reference_variant,
  active_records = nrow(target), active_keys_equal = keys_equal,
  max_probability_error = prob_error, probabilities_bitwise_equal = identical(r$probability, t$probability),
  max_pvalue_error = p_error, pvalues_exactly_equal = p_equal,
  significant_disagreements_lt005 = disagreement, significant_disagreements_le005 = disagreement_le,
  active_edge_counts_equal = edges_equal, graph_equal = graph_equal,
  graph_max_distance_error = graph_max_error, status = if (pass) 'PASS' else 'FAIL')
write.table(row, file.path(target_dir, 'validation.tsv'), sep = '\t', quote = TRUE, row.names = FALSE)
print(row)
if (!pass) stop('Numerical consistency failed: ', target_dir)
