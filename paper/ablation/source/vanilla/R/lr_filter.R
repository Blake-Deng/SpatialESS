#' Filter CellChat interactions by measured gene coverage
#'
#' Retains an interaction only when both ligand and receptor resolve to
#' non-empty gene sets and every required complex subunit occurs in `genes`.
#' This checks measured gene identifiers, not expression abundance, spatial
#' variability, overexpression or statistical significance. Zero-expression
#' rows are still measured genes. Cofactor genes are reported but are optional,
#' following the existing CellChat-compatible cofactor handling.
#'
#' @param lr CellChat interaction data frame, in the desired original order.
#' @param genes Unique, non-empty measured gene identifiers (normally symbols).
#' @param complex CellChat complex table with row names and subunit columns.
#' @param cofactor Optional cofactor table from the same database as `lr`.
#' @param ligand_col,receptor_col Names of the ligand and receptor columns.
#' @return A list with `lr` (retained rows), `excluded` (excluded original rows),
#'   `valid` (logical mask in input order), `audit` (one row per input LR), and
#'   `missing_genes` (one row per absent gene and component role). Missing-gene
#'   rows marked `required = FALSE` do not exclude an interaction. Empty or
#'   unresolved definitions are described in `audit$reason`. Original LR row
#'   names, columns and order are preserved. A warning reports all exclusions;
#'   the function does not silently discard records.
#' @export
#' @examples
#' lr <- data.frame(ligand = c("L", "L"), receptor = c("R", "RC"))
#' complex <- data.frame(subunit_1 = "R", subunit_2 = "R_missing",
#'                       row.names = "RC")
#' # The intentional warning reports the incomplete second interaction.
#' covered <- suppressWarnings(filter_cellchat_lr(lr, c("L", "R"), complex))
#' covered$valid
#' covered$missing_genes
filter_cellchat_lr <- function(lr, genes, complex, cofactor = NULL,
                              ligand_col = "ligand", receptor_col = "receptor") {
  report <- .spatialess_lr_audit(lr, genes, complex, cofactor,
                               ligand_col, receptor_col)
  .spatialess_warn_lr_exclusions(report)
  report
}

.spatialess_lr_audit <- function(lr, genes, complex, cofactor = NULL,
                                ligand_col = "ligand", receptor_col = "receptor") {
  lr <- as.data.frame(lr, stringsAsFactors = FALSE)
  genes <- as.character(genes)
  if (!length(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("genes must be unique, non-missing, non-empty identifiers.", call. = FALSE)
  }
  if (!all(c(ligand_col, receptor_col) %in% names(lr))) {
    stop("lr is missing ligand or receptor columns.", call. = FALSE)
  }
  subunits <- grep("^subunit", colnames(complex))
  if (is.null(rownames(complex)) || !length(subunits)) {
    stop("complex requires row names and subunit columns from the same database as lr.", call. = FALSE)
  }
  cof_cols <- if (!is.null(cofactor)) grep("cofactor", colnames(cofactor)) else integer()
  if (!is.null(cofactor) && (is.null(rownames(cofactor)) || !length(cof_cols))) {
    stop("cofactor requires row names and cofactor columns.", call. = FALSE)
  }
  resolve <- function(name, table, columns) {
    name <- as.character(name)
    if (length(name) != 1L || is.na(name) || !nzchar(name)) return(character())
    values <- if (name %in% rownames(table)) {
      unlist(table[name, columns, drop = FALSE], use.names = FALSE)
    } else name
    values <- as.character(values)
    values[!is.na(values) & nzchar(values)]
  }
  lig <- lapply(lr[[ligand_col]], resolve, table = complex, columns = subunits)
  rec <- lapply(lr[[receptor_col]], resolve, table = complex, columns = subunits)
  miss_l <- lapply(lig, setdiff, y = genes)
  miss_r <- lapply(rec, setdiff, y = genes)
  valid <- lengths(lig) > 0L & lengths(rec) > 0L &
    lengths(miss_l) == 0L & lengths(miss_r) == 0L
  interaction <- if ("interaction_name" %in% names(lr)) as.character(lr$interaction_name) else rownames(lr)
  join <- function(sets) vapply(sets, paste, character(1), collapse = ";")
  reasons <- vapply(seq_len(nrow(lr)), function(i) {
    if (valid[[i]]) return("complete")
    paste(c(if (!length(lig[[i]])) "unresolved_ligand_definition",
            if (!length(rec[[i]])) "unresolved_receptor_definition",
            if (length(miss_l[[i]])) "missing_ligand_genes",
            if (length(miss_r[[i]])) "missing_receptor_genes"), collapse = ";")
  }, character(1))
  audit <- data.frame(
    input_row = seq_len(nrow(lr)), interaction_name = interaction,
    ligand = as.character(lr[[ligand_col]]), receptor = as.character(lr[[receptor_col]]),
    ligand_genes = join(lig), receptor_genes = join(rec),
    missing_ligand_genes = join(miss_l), missing_receptor_genes = join(miss_r),
    valid = valid, reason = reasons, stringsAsFactors = FALSE
  )
  missing <- data.frame(input_row = integer(), interaction_name = character(),
                        role = character(), component = character(), gene = character(),
                        required = logical(), stringsAsFactors = FALSE)
  append_missing <- function(sets, names, role, required) {
    sizes <- lengths(sets)
    i <- rep(seq_along(sets), sizes)
    if (!length(i)) return(missing[FALSE, ])
    data.frame(input_row = i, interaction_name = interaction[i], role = role,
               component = as.character(names)[i], gene = unlist(sets, use.names = FALSE),
               required = required, stringsAsFactors = FALSE)
  }
  missing <- rbind(missing,
                   append_missing(miss_l, lr[[ligand_col]], "ligand", TRUE),
                   append_missing(miss_r, lr[[receptor_col]], "receptor", TRUE))
  if (!is.null(cofactor)) {
    for (role in intersect(c("co_A_receptor", "co_I_receptor", "agonist", "antagonist"), names(lr))) {
      sets <- lapply(lr[[role]], resolve, table = cofactor, columns = cof_cols)
      missing <- rbind(missing, append_missing(lapply(sets, setdiff, y = genes),
                                               lr[[role]], role, FALSE))
    }
  }
  list(lr = lr[valid, , drop = FALSE], excluded = lr[!valid, , drop = FALSE],
       valid = unname(valid), audit = audit, missing_genes = missing)
}

.spatialess_warn_lr_exclusions <- function(report) {
  n <- sum(!report$valid)
  if (n > 0L) warning(sprintf(
    "Excluded %d of %d LR records with incomplete ligand/receptor components; retained %d. See $audit and $missing_genes (or result$lr_filter when lr_missing = 'filter').",
    n, length(report$valid), sum(report$valid)), call. = FALSE)
  invisible(NULL)
}

.spatialess_lr_error <- function(report) {
  bad <- which(!report$valid)
  i <- bad[[1L]]
  detail <- report$missing_genes
  missing <- unique(detail$gene[detail$input_row == i & detail$required])
  stop(sprintf(paste0(
    "%d LR records have unresolved ligand/receptor components (first row: %d). ",
    "First interaction: %s; %s. Use filter_cellchat_lr(lr, rownames(expression), complex, cofactor), ",
    "then pass its $lr with the same database tables, or explicitly set lr_missing = 'filter'. ",
    "Check gene-symbol conventions; missing measured genes cannot be restored by sampling more cells."),
    length(bad), i, report$audit$interaction_name[[i]],
    if (length(missing)) paste("missing genes:", paste(missing, collapse = ", ")) else report$audit$reason[[i]]),
    call. = FALSE)
}
