# Panel coverage and the SpatialESS 0.1.2 interface

The default remains strict: `spatialess()` rejects an incomplete ligand or
receptor before inference. Use a coverage report explicitly:

```r
db <- CellChatDB
coverage <- filter_cellchat_lr(db$interaction, rownames(expression),
                              db$complex, db$cofactor)
result <- spatialess(expression, coordinates, group,
                     lr = coverage$lr, complex = db$complex,
                     cofactor = db$cofactor, nboot = 100, seed.use = 1)
```

Or deliberately request the same filtering inside the call:

```r
result <- spatialess(expression, coordinates, group,
                     lr = db$interaction, complex = db$complex,
                     cofactor = db$cofactor, nboot = 100, seed.use = 1,
                     lr_missing = "filter")
result$lr_filter$audit
result$lr_filter$missing_genes
```

Both paths emit a warning when exclusions occur. Neither silently edits the
original database. The audit preserves original input-row indices even when
the inference engine subsequently orders LR by signaling category. An LR
with an empty definition is rejected explicitly. If no LR survives, inference
stops with an actionable message; the standalone filter still returns its report.

## What coverage means

- Single ligand/receptor: its gene identifier must occur in the expression rows.
- Complex ligand/receptor: every non-empty listed subunit must occur.
- Zero-expression measured rows are retained. Missing rows are not imputed.
- Cofactor availability is recorded with `required = FALSE`. The established
  available-gene cofactor behavior is preserved; cofactors do not become new
  mandatory filters.
- This does not identify spatially variable, overexpressed or significant LR.
  Official V3's selected `chat@LR$LRsig` is another valid input.

Use one database object for all tables. Do not compute a mask on `CellChatDB`
and then apply it to a different `CellChatDB.human` object. A logical mask alone
does not filter anything: pass `coverage$lr` (or `lr[valid, , drop = FALSE]`).

For official comparisons preserve the entire prepared signaling expression
matrix, including its maximum used in normalization; do not trim the input
to only the filtered LR genes. Both engines must share this matrix and all
parameters, cells, groups and LR. Using different preparations can change
probabilities independently of the coverage interface.

## Real cervical dataset diagnosis

The 10x Xenium Prime FFPE Human Cervical Cancer H5 contains 840,387 cells,
9,475 features, and exactly 5,101 Gene Expression features. With the V3 database
(3,234 interactions), 1,076 LR have complete mandatory genes and 2,158 do not.
The first failing original row is TGFB1_ACVR1C_TGFBR2; receptor subunit ACVR1C
is absent from the full panel. Increasing the sampled cell count cannot restore
an unmeasured panel gene. This is distinct from the official sparse-distance
matrix dimension bug affecting trailing isolated cells.

The release includes a real 1k fixture under `validation/fixtures/`. It uses
23 official graph clusters, not manual cell-type annotations, with proportional
random sampling from cells satisfying transcript_counts >= 20 and valid
coordinates/cluster assignments. Gene counts were normalized to 10,000 per
cell and log1p transformed; all 5,101 measured rows were retained. The fixture
tests the interface and is not a new biological-validation claim.

Source: https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-cervical-cancer

## Regression and timing definitions

`validation/lr_interface/new_vs_old.tsv` compares the same input in isolated
R processes installed from old 0.1.1 and new 0.1.2 sources. The three modes are
old manual filtering, new explicit filtering, and new in-call filtering. The
manual coverage rule is independently frozen in the benchmark script; its
retained LR must equal the original prepared LR and each new selection exactly.
Tests compare input contracts, final RNG, complete result records, diagnostics,
probabilities, p-values, significant sets and result parameters.

The scale ladder uses the archived CosMx nested centre-region subsets at
1k/5k/10k/25k/50k/100k/200k/443,515, 373 LR and 20 permutations. A 10k
cluster-stratified CosMx subset is separate. Cervical 1k/5k tests use 1,076 LR
and 100 permutations; they do not reproduce the collaborator's private QC,
annotation or exact selected cells. Seeds and case-specific settings are saved.

`resources.tsv` records whole-process wall time/peak RSS (`/usr/bin/time -v`)
and timed inference. New-mode inference intervals include the added filter
and audit; old-mode intervals start after manual preparation. Whole-process
times include input loading and diagnostic output. These single runs verify
resource feasibility, not a precise speed comparison; shared-server workloads
and concurrent build/check activity were not controlled. The new audit has
measurable overhead. It does not imply the unchanged C++ core is slower.

No new official inference was launched for this interface-only release.
`new_vs_archived_official.tsv` compares the new output with existing official
V3 RDS results, after checking case size, LR count, permutations and spatial
parameters and reproducing the original prepared input and seed. Original and
dimension-patched comparators are separate labels. The original official
100k+ failures remain in the historical status tables and are not counted as
successful fidelity comparisons. The strict numerical contract uses RDS:
decimal TSV exports are for inspection and can round floating-point values.

Re-run from the repository root:

```bash
R CMD INSTALL .
Rscript benchmarks/verify_lr_interface_mini.R .
Rscript benchmarks/verify_gse250346_result_bundle.R results/tables
```

The scale runner accepts an installed library, case name (`cosmx`,
`cosmx_stratified`, `cervical`), cell count, mode and output directory. Set
`LR_REGRESSION_COSMX`, `LR_REGRESSION_STRATIFIED`, `LR_REGRESSION_CERVICAL` to
the respective downloaded/prepared data roots on another server.
