# Main Validation Record

This document separates exact numerical comparison from scalability and
biological interpretation.

## Evidence map

| Dataset | Role | Comparison scope | Status |
|---|---|---|---|
| GSE250346 | Primary multi-patient biology | 45 slices, 35 patients, 37 computable LR, patient-random-intercept LMM and slice-level GLM comparator | Complete with fit diagnostics |
| GSE306130 | Independent fidelity pilot | 4 samples, 1,000 cells/sample, 33 LR, 20 permutations | Complete pilot |
| GSE313006 | Independent fidelity fixture | 1,000 cells, 133 LR, 20 and 100 permutations | Complete fixture |
| CosMx Cancerous Liver | Paired scalability benchmark | 250, 1,000, 5,000, 10,000, 25,000 and 50,000 cells | Complete paired ladder |
| CosMx Cancerous Liver | Large-scale continuation | 100,000, 200,000 and 443,515 cells with SpatialESS | Reconciled: SpatialESS completed; official failed at all three scales |
| Xenium ovarian | Million-cell stress test | 1,156,091 cells, 648 genes, 2,357 groups, 627 LR | Complete engineering run; version audit required |

The compact table is results/tables/dataset_validation_overview.tsv.

## Numerical fidelity

The criteria are:

    active-record Jaccard = 1
    max absolute probability difference at floating-point scale
    exact p-value fraction = 1 for the same permutation design
    significance threshold disagreements = 0

The CosMx paired ladder satisfies these criteria at every paired size through
50,000 cells. GSE306130 and GSE313006 provide independent smaller fixtures.
GSE250346 provides a broad 45-sample check under a distinct annotation and
sample structure.

### GSE250346 official boundary

The original official run completed 41/45 samples. TILD117LA then completed
with the official implementation's minimal repair:

    Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(NC, NC))

This only preserves zero-degree cells that would otherwise be omitted by
automatic dimension inference. The repaired official output and SpatialESS
had 21,998 records, active Jaccard 1, maximum probability difference
7.81e-18, exact p-value fraction 1, and zero significance disagreements.

The remaining three samples are a separate official scale boundary:

| Sample | First limiting LR | Estimated nonzero outer product | Official outcome |
|---|---|---:|---|
| VUILD106MA | FN1-CD44 | 4,955,880,084 | signal 11 in Matrix::crossprod() |
| VUILD110LA | FN1-CD44 | 4,386,914,951 | signal 11 in Matrix::crossprod() |
| VUILD115MA | FN1-CD44 | 2,617,359,480 | signal 11 in Matrix::crossprod() |

Each exceeds the approximately 2.147 billion nonzero-index limit of the
32-bit dgCMatrix representation. This is an official implementation boundary,
not an input-data failure or a generic OOM claim.

## Runtime and memory

The paired CosMx ladder is the primary runtime/RSS figure source. At 50,000
cells, the recorded values were 23,919 seconds and 60.78 GiB for official
SpatialCellChat V3 versus 10.48 seconds and 0.802 GiB for SpatialESS. The official comparator did not complete at 100,000 cells and above. Failure modes are retained in results/tables/cosmx_paired_scaling.tsv: signal 11 at 100,000 cells, followed by explicit vector-allocation failures at 200,000 and 443,515 cells. SpatialESS completed all three larger fixtures with exit status 0. These scales are therefore SpatialESS scalability results, not paired numerical-fidelity comparisons.

The Xenium stress test is separate: it demonstrates that the SpatialESS engine
can operate on 1.156 million cells and 2,357 spatial groups. It should not be
used as a direct official comparison unless the current release commit and
parameters are re-audited against the historical source filenames.

## Multi-patient linear mixed model

The primary GSE250346 cohort analysis keeps one communication profile per spatial slice and uses patient as a random intercept:

    log1p(score) ~ condition + (1 | patient)

This model is fitted feature by feature with `lme4::lmer()`. It retains multiple slices from one patient while accounting for within-patient correlation. The simple slice-level model, `log1p(score) ~ condition`, is retained as a comparator. A patient-aggregated GLM remains an additional sensitivity analysis. Cells are never treated as independent cohort replicates.

The implementation standardizes each feature response before fitting and transforms estimates and standard errors back to the original log1p-score scale. It first uses `bobyqa`, retries numerical failures with `nloptwrap`, and records random-effect variance, residual variance, maximum gradient, optimizer, convergence messages and fit errors. Wald normal-approximation p-values are explicit; BH correction is restricted to `status = "ok"`.

Fit states are interpreted as follows:

- `ok` with `singular = FALSE`: converged model with nonzero patient variance;
- `ok` with `singular = TRUE`: converged boundary model with patient variance estimated as zero;
- `convergence_warning`: retained for diagnosis and excluded from FDR;
- `fit_failed`: both optimizers failed;
- `insufficient_nonzero`: not fitted because the feature occurs in too few slices.

Singular fits are not counted as software failures. The frozen run contained 39,242 `ok` fits (21,934 non-singular and 17,308 singular), 2,299 convergence warnings, 305 failures after both optimizers and 8,194 features below the nonzero-slice threshold. The LMM identified 5,158 FDR-significant features, compared with 4,172 for the slice-level GLM; 4,109 were shared (Jaccard 0.787). Valid-fit effect estimates had Spearman rho 0.9967 and 99.28% direction agreement. Formal result counts are frozen in `results/tables/gse250346_lmm_model_summary.tsv`, with status and diagnostic tables beside it.

The main biological summaries include collagen and matrix-associated programs such as COL1A1-CD44, COL1A2-CD44, FN1-CD44, and additional VEGF, CCL, SPP1 and MHC-II programs. These are biologically plausible disease-associated programs, not experimental validation.

## Robustness already available

The GSE250346 robustness outputs include the patient-random-intercept LMM, its slice-level GLM comparator, leave-one-patient-out fits, patient mean, median and cell-weighted aggregation, and min.cells.sr values 5, 10 and 20.

The central stability result is direction concordance: all primary significant
effects retained their direction in the leave-one-patient-out fits. Exact
retention and overlap values remain in the original TSV files.

## What is not claimed

- Official numerical agreement is not biological truth.
- A SpatialESS-only run is not an exact official comparison.
- A four-sample GLM pilot is not a patient-cohort disease conclusion.
- The Xenium stress test is not a leakage-correction validation.
- Official signal 11 is not described as a biological or input-data error.

