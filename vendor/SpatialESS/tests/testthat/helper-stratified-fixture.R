make_stratified_analytic_fixture <- function(seed = 29L, nested = FALSE,
                                             tied = FALSE, mild_zero = FALSE) {
  set.seed(seed)
  block_count <- 8L
  cells_per_block <- 40L
  cell_count <- block_count * cells_per_block
  cell_id <- paste0("cell", seq_len(cell_count))
  block <- rep(seq_len(block_count), each = cells_per_block)
  position <- rep(seq_len(cells_per_block), block_count)
  group <- if (nested) {
    factor(ifelse(block <= block_count / 2L, "A", "B"),
           levels = c("A", "B"))
  } else {
    factor(ifelse(position %% 2L == 0L, "A", "B"),
           levels = c("A", "B"))
  }
  block_baseline <- stats::rgamma(block_count, shape = 3, rate = 2)
  ligand <- block_baseline[block] *
    stats::rgamma(cell_count, shape = 6, rate = 3)
  receptor <- block_baseline[block] *
    stats::rgamma(cell_count, shape = 7, rate = 3)
  if (tied) {
    ligand[position %% 4L != 0L] <- 0
    receptor[position %% 3L != 0L] <- 0
  }
  if (mild_zero) {
    zero_index <- position <= 8L
    ligand[zero_index] <- 0
    receptor[zero_index] <- 0
  }
  expression <- Matrix::Matrix(
    rbind(L = ligand, R = receptor), sparse = TRUE,
    dimnames = list(c("L", "R"), cell_id)
  )
  coords <- cbind(x = position, y = block * 100)
  rownames(coords) <- cell_id
  graph <- build_radius_graph_csr(
    coords, radius = 45, sample_id = rep("s1", cell_count),
    weight = "binary"
  )
  support <- build_group_support_csr(graph, group)
  blocks <- build_spatial_blocks(
    coords, block_size = c(100, 100), origin = c(0, 100)
  )
  complex <- data.frame(subunit_1 = character(), row.names = character())
  cofactor <- data.frame(cofactor1 = character(), row.names = character())
  lr <- data.frame(interaction_name = "L_R", ligand = "L", receptor = "R")
  components <- prepare_cellchat_lr_components(
    lr, rownames(expression), complex, cofactor
  )
  prepared <- prepare_sparse_trimean(
    expression, genes = components$genes, normalize = FALSE
  )
  list(group = group, support = support, blocks = blocks,
       components = components, prepared = prepared)
}
