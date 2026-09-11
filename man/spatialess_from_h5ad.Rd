\name{spatialess_from_h5ad}
\alias{spatialess_from_h5ad}
\title{Run SpatialESS from a SPARKLE-corrected H5AD file}
\usage{
spatialess_from_h5ad(
  h5ad, group_col, coordinate_cols, lr, complex, cofactor,
  gene_list = NULL, normalization = c("log1p_cp10k", "none"),
  metadata_file = NULL, cell_col = "cell", layer = NULL,
  python = Sys.which("python3"), bridge_script = NULL,
  output_dir = NULL, keep_bundle = FALSE, filter_unresolved = TRUE,
  lr_missing = NULL, ...
)
}
\arguments{
\item{h5ad}{Path to a SPARKLE-corrected H5AD file.}
\item{group_col}{Metadata column containing cell-group labels.}
\item{coordinate_cols}{Two or three spatial-coordinate column names.}
\item{lr}{CellChat-compatible ligand-receptor interaction table.}
\item{complex}{CellChat-compatible complex table.}
\item{cofactor}{CellChat-compatible cofactor table.}
\item{gene_list}{Optional genes to export; by default extracted from the database.}
\item{normalization}{Either \code{"log1p_cp10k"} or \code{"none"}.}
\item{metadata_file}{Optional TSV containing labels and coordinates.}
\item{cell_col}{Cell identifier column in \code{metadata_file}.}
\item{layer}{Optional H5AD layer; \code{NULL} uses \code{X}.}
\item{python}{Python executable containing the required dependencies.}
\item{bridge_script}{Optional explicit path to the Python bridge script.}
\item{output_dir}{Optional persistent bridge-bundle directory.}
\item{keep_bundle}{Whether to retain an automatically created bundle.}
\item{filter_unresolved}{Legacy policy: TRUE filters incomplete LR with a warning and full audit; FALSE uses strict checking.}
\item{lr_missing}{Explicitly select \code{"filter"} or \code{"error"}. NULL preserves \code{filter_unresolved}; conflicting explicit options are rejected.}
\item{...}{Additional arguments passed to \code{SpatialESS::spatialess()}.}
}
\description{
Converts a corrected H5AD file without densifying it in R, removes LR rows
that cannot be resolved from the exported genes, and runs the exact
SpatialESS engine preserving the tested SpatialCellChat V3 definitions.
}
\value{A SpatialESSResult with bridge diagnostics. Filtering stores the complete input LR audit, excluded rows and missing genes in \code{result$lr_filter}. A persistent bundle also contains \code{lr_audit.tsv} and \code{lr_missing_genes.tsv}.}
