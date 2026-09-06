# SpatialESS 0.1.0 Release Record

Release date: 2026-09-06

This release is the reproducible software and benchmark snapshot for the
SpatialESS main method line. It is intended for independent verification of
the implementation, numerical fidelity and resource scaling. Figures and
large source datasets are intentionally outside this release snapshot.

## Frozen scope

- spatial_v3_compatible() is the official-compatible numerical-fidelity path.
- spatialess() is the integrated sparse spatial workflow.
- Main inference uses fixed permutation inference; experimental analytic and
  hybrid modes are retained but are not part of the primary claim.
- Exact paired comparisons use identical nested cells, coordinates, gene panel,
  group labels, LR definitions, seed and one-thread settings.

## Frozen parameters

| Setting | Value |
|---|---|
| Contact radius | 15 um effective radius |
| Diffusion radius | 35 um effective radius |
| Spatial block size | 100 um |
| Permutations | 20 for the compact benchmark tables |
| Native threads | 1 |
| BLAS/OpenMP threads | 1 |
| Official memory envelope in CosMx ladder | 120 GiB for official runs |
| Main significance threshold | 0.05 |

## Final status boundaries

- GSE250346 contains 45 samples from 35 patients: 9 controls and 26 disease
  samples. SpatialESS completed 45/45.
- The unmodified official GSE250346 workflow completed 41/45. TILD117LA
  completed after the minimal official dims = c(NC, NC) repair and passed
  exact comparison. Three larger samples reached the official sparse-matrix
  boundary at FN1-CD44 and exited with signal 11.
- CosMx paired fidelity is available from 250 through 50,000 cells. SpatialESS
  also completed 100,000, 200,000 and 443,515 cells. The official process
  failed at those three scales, so they are reported as SpatialESS-only
  scalability results.

See results/tables/cosmx_paired_scaling.tsv,
results/tables/cosmx_queue_reconciliation.tsv and
results/tables/gse250346_official_boundaries.tsv for the machine-readable
status records.

## Verification

From the repository root:

  /home/dzf/miniforge3/envs/cellchat-acceleration/bin/R CMD check --no-manual --no-build-vignettes SpatialESS_0.1.0.tar.gz

The expected package-check result for this snapshot is Status: OK.
Independent users should run the small numerical tests first and only then
provide their own prepared data roots for the larger benchmark drivers. Raw
data and benchmark RDS files are excluded from Git by design.
