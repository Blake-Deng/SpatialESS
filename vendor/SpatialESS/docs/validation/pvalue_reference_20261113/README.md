# SpatialESS 0.1.4.9003 reference-parity validation

The validated implementation preserves the strict `Pboot > Pobs`, `count/nboot`
empirical rule. No epsilon, p-value correction or replacement is used.

## Fresh official paired comparisons

Original SpatialCellChat 0.1.0 was rerun in the same server environment, pinned
to source commit `be88f300464cf970b36a6165ea7d2e41a47a7511`. Seven installed
core functions were compared with that source and had identical formals/bodies.

Six real cases and five synthetic cases produced **17,880 active records**.
The actual SpatialESS API's record domain, probabilities and p-values matched
the fresh official execution element by element using exact equality. Captured
actual permutation indices and permuted labels also matched.

- Lung VUILD142MA: 9,222 cells, 37 LR interactions, 20 permutations; 3,386 records.
- Lung TILD049MA: 11,140 cells, 37 LR interactions, 20 permutations; 10,924 records.
- Four GSE306130 cases: 1,000 cells each, 33 LR interactions, 20 permutations.
- Synthetic cases: single-gene ties (100 permutations), complexes/modifiers,
  percentage boundaries, mixed contact modes, and non-unit coordinate ratio.

See [exact_summary.tsv](exact_summary.tsv). The diagnosed VUILD142MA
COL1A2_CD44 / Monocytes/MDMs-to-Monocytes/MDMs p-value is 0.70 in the reference
and corrected implementation, compared with 0.65 in the previous engine.

For each captured LR/permutation, the production aggregation helper separately
replayed official cell probabilities and matched the full observed/null/count/p
arrays. This is helper replay, not direct interception of every internal C++
intermediate. Public tail counts, when absent, were reconstructed from p*nboot;
that reconstruction is an auxiliary diagnostic rather than independent evidence.

## Environment and checks

The paired server environment used R 4.5.3, Matrix 1.7-5, OpenBLAS 0.3.33 and
one BLAS/OMP thread. Both methods used identical input order, groups, parameters,
RNG settings and seeds. Exact cross-version or cross-BLAS equality is not asserted.

The final core passed 1,681 unit assertions locally, with zero failures, errors,
warnings or skips; see [unit_totals.tsv](unit_totals.tsv). Both server package
checks ended with Status: OK. The engine's server package check skipped tests,
which were run separately locally; the server bridge check ran its own tests.
The installed bridge smoke check passed 16 checks.

The algorithm/API/test core is unchanged by publishing documentation and workflow
updates. Canonical core-file hash: `bd9a0076a02e592f3097b6908ae5ac307b7a4cba47557687a4ca3c8b1c0da827`. It is SHA256 of a compact,
key-sorted JSON object mapping relative files under R/src/man/tests to each
file's SHA256, excluding compiled .o/.so files. Main and the bundled engine match.

## Historical audit and limits

A preceding candidate was also checked against 66 saved historical references
(57 unmodified official and 9 dimension-only repaired references). All input and
reference file MD5 values were checked. Twenty-two comparisons were fully exact;
40 lung comparisons had identical p-values but probability differences up to
8.67e-19, whose historical cause was not diagnosed entry by entry. Four small
cases were not strictly paired with the current group input: the prior engine
also disagreed with their historical records, while fresh official runs on the
current inputs matched the corrected engine exactly. TILD049MA, which had a
historical last-bit difference, likewise matched its fresh official run exactly.
This is not a claim of 66 fresh exact official reruns.

Reference mode requires dense cell-by-group/group-by-group workspaces. The
per-matrix entry guard is not a total-RAM bound. Historical stream-mode speed
and memory figures must not be relabeled as new reference-mode benchmarks.
Spatial-block, adaptive and hybrid inference remain separate experimental nulls.
No new biological-validation claim or DOI is made by this pre-release.

Historical paper/source directories retain research-wrapper metadata and shared
source aliases; they are not the current package installation entry point.
Install from the repository root (or SPARKLE vendor/SpatialESS). Reproduce old
versions from their corresponding archived tags and materials.
