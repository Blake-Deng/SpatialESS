# LR interface and paired main / SPARKLE repositories

The two repositories share the same SpatialESS 0.1.2 inference and GLM/LMM
implementation. Main takes an R matrix and spatial metadata. This integration
adds SPARKLE installation instructions and an H5AD-to-R conversion step,
then calls the bundled main engine. The upstream SPARKLE correction code and
Python H5AD conversion code have not been changed in this update.

Use `SpatialESS::filter_cellchat_lr()` for explicit coverage filtering, or
`lr_missing = "filter"` in either inference entry point. All mandatory complex
subunits must be present. Missing cofactors are recorded without changing their
established handling. Main defaults to `"error"`; the bridge preserves its
previous filtering default. The bridge's `filter_unresolved` argument remains
supported. Exclusions now generate a warning and a full input-row audit rather
than being silently swallowed by a per-LR exception handler. Other input errors
are allowed to fail explicitly.

The bridge returns `result$lr_filter` with `lr`, `excluded`, `valid`, `audit`
and `missing_genes`. If `output_dir` is supplied, it also writes `lr_audit.tsv`
and `lr_missing_genes.tsv`. These files record coverage of the exported matrix,
so a deliberately restricted `gene_list` can also cause exclusions. Do not
trim the normalization matrix to only retained LR genes when comparing engines.

## Validation design

The September 10 regression uses the previous ovarian RAW and corrected H5AD
inputs with shared RAW-derived RCTD labels and unchanged coordinates. The nested
cell subsets are the first 1,000, 5,000 and all 16,247 cells in archived order.
They are software regression inputs, not independent biological replicates.
Each condition and size is run in three isolated processes:

1. Old 0.1.1 bridge and 0.1.1 engine, with the existing coverage selection.
2. Updated bridge and bundled 0.1.2 engine, with explicit audited filtering.
3. Direct main 0.1.2 inference on the same bridge-exported R matrix, with
   explicit coverage filtering.

The frozen ovarian database has 1,939 LR candidates, of which 1,828 are
computable from 967 exported signaling genes. It is the database used by the
archived official **SpatialCellChat V3** run, loaded there from `CellChat`.
This is distinct from the newer 3,234-candidate database distributed with the
installed SpatialCellChat package. Algorithm generation and database release
are separate versioned inputs. Neither table is silently substituted in a
comparison. The new default CLI loads the SpatialCellChat database; reproducing
the old study requires the frozen database and gene list included with the mini.

All runs use 100 permutations, seed 20260903, one native/BLAS thread, ratio 1,
tol 5, interaction.range 30, contact.range 10, scale.distance 20,
min.percent 0 and min.cells.sr 1. Significance is `pvalue < 0.05`, as in the
archived ovarian study. The separate main regression records its own threshold
convention. No p-value definition is changed by this release.

`validation/bridge_interface/comparisons.tsv` reports exact records, input,
result-contract and RNG equality; probability differences, p-value agreement,
active Jaccard and significance disagreements. `resources.tsv` records wall
time, peak RSS and completion. Single timings on a shared server, including
occasional overlapping local checks, are not a controlled speedup experiment.
Direct-main timing excludes H5AD conversion; bridge timing includes it.

`history.tsv` separately compares full outputs with archived SpatialESS and
official V3 RDS results. Official inference and SPARKLE correction are not
rerun for this interface change. Upstream H5AD sources and archived outputs
are checksummed. Small floating-point differences from legacy text export are
reported rather than described as bitwise equality.

## Independent mini reproduction

From the repository root after installing the dependencies:

```sh
Rscript scripts/install.R .
PYTHON="$(command -v python)" Rscript benchmark/verify_ovarian_mini.R .
Rscript vendor/SpatialESS/benchmarks/verify_gse250346_result_bundle.R vendor/SpatialESS/results/tables
```

The real 1k RAW and corrected H5AD fixtures, frozen database, labels, gene list
and full-precision expected records are included under `validation/fixtures/`.
This mini verifies integration output; it does not refit SPARKLE from masks.
The original synthetic 24-cell demo is also retained.

The larger tests can be rerun with `benchmark/run_bridge_regression.R` and
the original study data. The script documents its arguments and data paths;
set `OVARIAN_SOURCE` for another machine. Full input matrices are not included
in the GitHub archive. Main's separate real cervical mini and 443,515-cell
regression evidence are included under `vendor/SpatialESS/validation/`.

SPARKLE is an independent upstream method. Its inclusion does not mean
calibration has been proven more biologically accurate in every dataset.
