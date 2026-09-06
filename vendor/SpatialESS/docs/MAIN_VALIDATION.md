# Main Validation Record

This document separates exact numerical comparison from scalability and
biological interpretation.

## Evidence map

| Dataset | Role | Comparison scope | Status |
|---|---|---|---|
| GSE250346 | Primary patient-level biology | 45 samples, 35 patients, 37 computable LR, sample-level records and patient GLM | Complete |
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

## Patient-level GLM

The primary GSE250346 analysis uses the patient, not the cell, as the unit of
inference:

    sample communication records
        -> patient-level aggregation
        -> condition and covariate model
        -> FDR correction across communication features
        -> conserved / disease-associated programs

The prepared cohort contains 45 samples, 35 patients, 9 controls and 26
disease patients. Covariate-aware and aggregation sensitivity outputs are
kept under the benchmark result root and summarized in
results/tables/patient_glm_design.tsv and
results/tables/patient_glm_top_programs.tsv.

The main biological summaries include collagen and matrix-associated programs
such as COL1A1-CD44, COL1A2-CD44, FN1-CD44, and additional VEGF, CCL, SPP1 and
MHC-II programs. These are biologically plausible disease-associated programs,
not experimental validation.

## Robustness already available

The GSE250346 robustness outputs include leave-one-patient-out fits, patient
mean, median and cell-weighted aggregation, min.cells.sr values 5, 10 and 20,
and mixed-model sensitivity fits with singular-fit diagnostics.

The central stability result is direction concordance: all primary significant
effects retained their direction in the leave-one-patient-out fits. Exact
retention and overlap values remain in the original TSV files.

## What is not claimed

- Official numerical agreement is not biological truth.
- A SpatialESS-only run is not an exact official comparison.
- A four-sample GLM pilot is not a patient-cohort disease conclusion.
- The Xenium stress test is not a leakage-correction validation.
- Official signal 11 is not described as a biological or input-data error.

