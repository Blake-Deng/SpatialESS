# Cervical Collaborator Diagnostic, 2026-09-11

## Important Result: Default-Filter Mismatch

Both official runs completed with NONZERO group-level results. However,
SpatialESS 0.1.2 does not exactly reproduce the official `min.percent=0.1`
filter at the rounding boundary in this submitted barcode subset:

| Shared input | LR | Official active | ESS active | Official p<0.05 | ESS p<0.05 | Significance disagreements |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Normalized, no ALRA/MERINGUE | 1076 | 57 | 56 | 31 | 32 | 3 |
| ALRA + MERINGUE | 596 | 232 | 230 | 123 | 125 | 2 |

Maximum probability differences on matched active records are 1.39e-17 and
2.08e-17 respectively. This is NOT a strict-fidelity pass: support and p-values
differ, despite the small probability differences on common records.

Concrete example: for `GJA1_GJA1`, Cluster_06 has 8/82 positive cells, or
0.09756098. Official `format(dataLR_temp[, -1], digits=1)` produces the string
`"0.10"`, which passes its subsequent comparison to 0.1. ESS compares the
unrounded fraction numerically and excludes the record. The same filtering
distinction can affect permutations and therefore other records' p-values.

The additional `isolate_percent_filter.R` script performs causal diagnostics:

- `numeric_normalized` and `numeric_alra_meringue` replace only this formatting
  step in a LOCAL COPY of the official group helper. These alter official
  filter semantics and must NOT be presented as fidelity to original V3.
- `relaxed_normalized` uses the unchanged dimension-patched official function
  with min.percent=0, min.cells.sr=1 in both engines, matching the earlier
  ladder's filtering convention.

Production packages and release ZIPs are unchanged. Before claiming default
V3 equivalence for this case, ESS needs an explicit compatibility decision and
regression coverage for the official percentage-formatting behavior.

All three causal checks completed:

| Check | Reference active | ESS active | Active Jaccard | Exact p fraction | Significance disagreements |
| --- | ---: | ---: | ---: | ---: | ---: |
| Numeric-filter diagnostic, normalized | 56 | 56 | 1 | 1 | 0 |
| Numeric-filter diagnostic, ALRA + MERINGUE | 230 | 230 | 1 | 1 | 0 |
| Unchanged V3 helper, shared relaxed thresholds | 364 | 364 | 1 | 1 | 0 |

For the relaxed-threshold check, max absolute probability difference is
2.775558e-17 and both engines have 150 records with p<0.05. This confirms the
earlier benchmark convention on these new barcodes, but does not cure the
default-threshold mismatch. In the two numeric-filter diagnostic checks,
changing only the percentage-formatting rule removes ALL observed support,
p-value and significance differences in these reconstructed inputs.

## Scope

This is a controlled reconstruction, not a claim that the collaborator's full
session has been reproduced. The submitted `benchmark_cells.rds` contains only
barcode vectors (`n1000`, `n5000`, `n10000`, `n20000`, `n50000`). It does not contain
expression, coordinates, cell-type annotation, the official object, or results.

The 1,000 submitted barcodes and their order are preserved exactly. Expression
and coordinates are reconstructed from the already downloaded 10x Xenium Prime
FFPE Human Cervical Cancer data. Public 10x graph clusters (23 groups) substitute
for the unavailable collaborator cell-type labels. They are NOT biological
annotations. Counts are normalized by log1p(CP10K) over all 5,101 measured genes.

Source: https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-cervical-cancer

## Confirmed Observations

- Every submitted 1k barcode matches the public data. No barcode is duplicated.
- At ratio=1 and interaction.range+tol=35 um there are 46 directed non-self
  edges (23 unordered pairs). 954 sampled cells have no other sampled neighbor.
- The inferred original distance matrix is 927 x 927; specifying
  `dims=c(NC,NC)` restores 1000 x 1000.
- The minimum positive distance is 5.449088912 um. At scale.distance=0.01,
  the official distance check rejects 0.054490889 < 1. This is an ERROR, not
  a valid all-zero result. The exact error is saved in
  `results/default_scale_official_error.txt`.
- Official `preProcessing()` performs ALRA imputation and overwrites
  `data.signaling`. Its input must be shared with ESS if it is used.
- `identifyOverExpressedInteractions()` populates `LR$LRsig`. The official
  inference function uses that table unless `LR.use` is explicitly supplied.
- `computeCommunProb()` returns cell-level `net$prob.cell` and caches, not
  group-level `net$prob` and `net$pval`. Run `computeAvgCommunProb()` before
  counting group-level significant records.
- R's `$` can partially match `prob` to `prob.cell`. Use exact `[["prob"]]`
  and `[["pval"]]`, and fail if either is NULL. `sum(NULL < 0.05)` is zero,
  which is a missing-stage artifact, not evidence that probabilities are zero.

These facts do NOT prove which issue caused the collaborator's reported zero.
The original `data_1k`, `chat_1k`, ESS result and downstream commands remain
necessary for that conclusion. The reconstructed MERINGUE run selects 199
features and 596 LR, not their reported 190 and 530; do not conflate them.

## Controlled Comparisons

Two official runs use the same matrix, barcode order, public group labels,
LR list and parameters as their corresponding ESS run:

1. `normalized`: no ALRA or MERINGUE; all 1076 covered LR.
2. `alra_meringue`: ALRA then MERINGUE; the resulting 596 LR in both engines.

Common parameters: ratio=1, tol=5, interaction.range=30, contact.range=10,
scale.distance=1, Kh=0.5, n=1, use.AGAN=TRUE, contact.dependent=TRUE,
contact.dependent.forced=FALSE, min.percent=0.1, min.cells.sr=5, avg.type=avg,
nboot=100, seed.use=1, colocalization.use=FALSE.

Scale.distance=1 passes the official lower-bound check. It is a shared diagnostic
setting, not a claim to reproduce the collaborator's unknown parameter setting.
Official cell-level work uses 8 multisession workers; group/permutation stages
are sequential. Native thread counts are 1. These are diagnostics, NOT an
isolated runtime benchmark.

SpatialESS 0.1.2 is loaded from the verified release library without source edits.
Official V3 is the independent dimension-patched copy from commit
`be88f300464cf970b36a6165ea7d2e41a47a7511`. The only official source change is
the explicit distance-matrix dimension. All raw inputs remain unmodified.

Each `results/official_<condition>/` contains:

- `chat_cell_stage.rds`, `cell_stage.tsv`: state immediately after cell inference.
- `chat_complete.rds`: state after averaging and permutation.
- `records_comparison.tsv`, `comparison.tsv`: keyed numerical comparison.
- `provenance.txt`, `sessionInfo.txt`, `COMPLETED.txt`: environment and status.

The authoritative state is the presence of a completed comparison, not this
README alone. Do not label pending or failed runs as passed.

## Next Input From Collaborator

Source `export_collaborator_debug.R` in their R session, then run:

```r
export_collaborator_debug(
  data_1k, chat_1k, result_1k_ess,
  out = "cervical_debug"
)
```

It writes a NEW directory and refuses to overwrite an existing directory.
Also supply the actual `computeCommunProb1()` invocation, the subsequent
aggregation commands, and the exact command used to obtain zero.
Defining `computeCommunProb1` does not replace `computeCommunProb` in a package
namespace; verify which function was actually called.

## Reproduction

Use `/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript` and the existing
package libraries recorded by the scripts. First run `export_cells.R`, then
`python prepare_counts.py`, then `prepare_and_diagnose.R`. Run
`run_official.R normalized` and `run_official.R alra_meringue` after preparation.
`check_default_guard.R` captures the official scale.distance rejection.

These scripts intentionally reproduce inputs rather than modify either engine.
