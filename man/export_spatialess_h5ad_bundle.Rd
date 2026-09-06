\name{export_spatialess_h5ad_bundle}
\alias{export_spatialess_h5ad_bundle}
\title{Export an H5AD expression matrix for SpatialESS}
\usage{
export_spatialess_h5ad_bundle(
  h5ad, group_col, coordinate_cols, gene_list = NULL, output_dir,
  normalization = c("log1p_cp10k", "none"), metadata_file = NULL,
  cell_col = "cell", layer = NULL, python = Sys.which("python3"),
  bridge_script = NULL, overwrite = FALSE
)
}
\arguments{
\item{h5ad}{Path to a SPARKLE-corrected H5AD file.}
\item{group_col}{Metadata column containing cell-group labels.}
\item{coordinate_cols}{Two or three spatial-coordinate column names.}
\item{gene_list}{Optional character vector of genes to export.}
\item{output_dir}{Destination directory for the sparse bridge bundle.}
\item{normalization}{Either \code{"log1p_cp10k"} or \code{"none"}.}
\item{metadata_file}{Optional TSV containing labels and coordinates.}
\item{cell_col}{Cell identifier column in \code{metadata_file}.}
\item{layer}{Optional H5AD layer; \code{NULL} uses \code{X}.}
\item{python}{Python executable containing the required dependencies.}
\item{bridge_script}{Optional explicit path to the Python bridge script.}
\item{overwrite}{Whether to replace a non-empty output directory.}
}
\description{
Exports SPARKLE-corrected expression, shared labels and coordinates as a
validated sparse SpatialESS bundle. CP10K denominators are computed over all
H5AD genes before optional signaling-gene selection.
}
\value{The normalized bundle path, invisibly.}
