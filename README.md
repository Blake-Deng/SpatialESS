# SpatialESS + SPARKLE

SpatialESSSPARKLE 0.1.4 connects sparse H5AD expression to SpatialESS 0.1.4. The complete engine, including patient-aware GLM/LMM analysis, is bundled in `vendor/SpatialESS`.

SPARKLE is independent upstream RNA-leakage correction software. This repository supplies the bridge and inference workflow; see [attribution](docs/SPARKLE_ATTRIBUTION.md).

## Method boundary

```text
Raw spatial data + segmentation/out-of-mask evidence
                 |
                 v
       SPARKLE (external software)
                 |
                 v
     corrected sparse H5AD expression
                 |
                 v
 SpatialESSSPARKLE H5AD bridge
  - align cells, labels and coordinates
  - compute all-gene CP10K denominators
  - export only required signaling genes
                 |
                 v
 SpatialESS inference
  - CSR radius graph
  - streamed LR and permutation computation
                 |
                 v
 sender-receiver-LR probabilities and p-values
```

The ovarian scFFPE reference and its official barcode annotation used in the
validation study are tissue-matched benchmark resources. They are not required
for users whose H5AD already contains cell-type labels.

## Installation

A C++17 compiler, R >= 4.1 and Python >= 3.9 are required. The validated SpatialESS engine is bundled in this repository and installed locally by the setup script.

```bash
git clone https://github.com/Blake-Deng/SpatialESS-SPARKLE.git
cd SpatialESS-SPARKLE

conda env create -f environment.yml
conda activate spatialess-sparkle
Rscript scripts/install.R .
```

The environment provides the Python dependencies and the pinned upstream
SPARKLE package, which is run before this bridge. The minimal example needs
neither CellChat nor SpatialCellChat. For the default command-line database,
install [SpatialCellChat](https://github.com/jinworks/SpatialCellChat).
The R API also accepts existing compatible database tables. Its legacy
`gene_list = NULL` path uses CellChat solely to extract database gene names;
passing an explicit `gene_list` avoids that optional dependency. Database
versions and normalization gene sets must be frozen for numerical comparisons.

## Minimal reproducible example

```bash
PYTHON="$(command -v python)" \
RSCRIPT="$(command -v Rscript)" \
bash examples/minimal/run_demo.sh
```

The example creates a sparse SPARKLE-like H5AD, converts it, runs one synthetic
LR pair with ten deterministic permutations, and writes:

```text
examples/minimal/output/demo_records.tsv
examples/minimal/output/demo_result.rds
```

## One-command analysis

When `cell_type`, `x_um` and `y_um` are already present in H5AD `obs`:

```bash
Rscript scripts/run_spatialess_from_h5ad.R \
  corrected.h5ad cell_type x_um y_um spatialess_result.rds
```

The command uses `CellChatDB.human`, 100 permutations and seed 1. For explicit
control over the database and SpatialESS parameters, use the R API:

If R and Python are in different environments, prefix the command with
`PYTHON=/absolute/path/to/python`; that interpreter must have the bridge
dependencies installed.

```r
library(SpatialCellChat)
library(SpatialESSSPARKLE)

data(CellChatDB.human, package = "SpatialCellChat")

result <- spatialess_from_h5ad(
  h5ad = "corrected.h5ad",
  group_col = "cell_type",
  coordinate_cols = c("x_um", "y_um"),
  lr = CellChatDB.human$interaction,
  complex = CellChatDB.human$complex,
  cofactor = CellChatDB.human$cofactor,
  gene_list = SpatialCellChat::extractGene(CellChatDB.human),
  lr_missing = "filter",
  normalization = "log1p_cp10k",
  nboot = 100L,
  seed.use = 1L
)

head(result$records)
result$bridge
result$lr_filter$audit
result$lr_filter$missing_genes
```

Incomplete ligand/receptor complexes require **all** mandatory subunits.
Excluded LR produce a warning and remain listed in the audit, with original
input-row indices and missing genes. With `output_dir` set, the bridge also
writes `lr_audit.tsv` and `lr_missing_genes.tsv`. Cofactor absences are reported
separately and keep the established optional-gene behavior. Use
`lr_missing = "error"` for strict checking. The bridge filters unresolved LR by default; the standalone main API retains its strict default.
The older `filter_unresolved` option remains supported.

For main-versus-bridge comparisons, use the **same exported matrix** from
`read_spatialess_bundle()`, LR, labels, coordinates, parameters and seed. Do
not remove extra signaling rows after export: they can set the normalization
maximum in inference. RAW and corrected inputs can produce different results;
that is separate from software numerical fidelity.

## External annotation

If labels or coordinates are not stored in H5AD, provide a tab-separated file:

```text
cell    shared_cell_type    x_um    y_um
cell_1  Fibroblast          120.3   88.1
cell_2  Macrophage          126.8   91.4
```

```r
result <- spatialess_from_h5ad(
  h5ad = "corrected.h5ad",
  metadata_file = "annotation.tsv",
  cell_col = "cell",
  group_col = "shared_cell_type",
  coordinate_cols = c("x_um", "y_um"),
  lr = CellChatDB.human$interaction,
  complex = CellChatDB.human$complex,
  cofactor = CellChatDB.human$cofactor,
  gene_list = SpatialCellChat::extractGene(CellChatDB.human),
  lr_missing = "filter",
  nboot = 100L,
  seed.use = 1L
)
```

Every H5AD cell must occur exactly once in the metadata. The bridge restores the
H5AD cell order and fails on missing labels, duplicate identifiers, non-finite
coordinates or unresolved matrix dimensions.

## Input assumptions

- Use `normalization = "log1p_cp10k"` only when the selected H5AD matrix or
  layer contains non-negative count-scale values.
- Use `normalization = "none"` when the matrix has already received the intended
  normalization.
- Coordinate units must match `interaction.range`, `contact.range`, `tol` and
  `scale.distance`.
- Cell-type annotation must be biologically appropriate and consistent across
  samples; this package does not infer it.
- SPARKLE requires its own platform-specific segmentation and out-of-mask
  evidence. Follow the upstream SPARKLE documentation before running this
  bridge.


## Cohort analysis and results

Use the bundled `SpatialESS::prepare_multisample_communication()` and `SpatialESS::fit_multisample_communication_lmm()` for multi-patient analysis. See [model documentation](vendor/SpatialESS/docs/MULTISAMPLE_LMM.md).

Current RAW and corrected-input bridge measurements are in `results/tables/bridge_results.tsv`. Direct-inference records, official comparator metrics and cohort fits are indexed by the [bundled result catalog](vendor/SpatialESS/docs/RESULTS.md).

Local checks are in `results/checks.tsv`. GitHub cloud Actions are **NOT RUN** for this source ZIP.
