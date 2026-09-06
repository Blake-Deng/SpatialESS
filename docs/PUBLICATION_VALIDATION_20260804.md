# Publication validation protocol

Date: 2026-08-04

## Scope

The formal comparator is SpatialCellChat V3. Numerical fidelity is evaluated
only for the `spatial_v3_compatible()` path. The integrated `spatialess()`
workflow is evaluated for deterministic execution, resource stability,
permutation resolution, radius sensitivity, real-data application and scaling.

No independent spatial CCC method is treated as numerical ground truth.

## Fixed environment

- one benchmark process at a time;
- one native and one BLAS thread;
- the same nested CosMx CancerousLiver cells and 373 computable LR;
- binary CSR graph weights;
- 15 um contact radius, 35 um diffusion radius and 100 um blocks unless varied;
- finite-corrected, greater-than-or-equal empirical p-values;
- every significant result fixed-permutation confirmed.

## Additional experiments

1. **Reproducibility:** five independent 5,000-cell processes with identical
   parameters and seed. Report exact probability/p-value fractions, threshold
   disagreements, runtime coefficient of variation and peak-RSS coefficient of
   variation.
2. **Permutation sensitivity:** 20, 50, 100 and 200 permutations at 5,000
   cells. Report runtime/RSS, p-value agreement with the 200-permutation
   reference, and significance-set Jaccard. Probabilities must remain exact
   because permutation count changes only null estimation.
3. **Seed stability:** five independent 100-permutation runs with different
   seeds. Report p-value rank concordance and significance-set overlap.
4. **Radius sensitivity:** diffusion radii 20, 35, 50, 75 and 100 um with the
   contact radius fixed at 15 um. Report graph edges, runtime/RSS, probability
   rank concordance and significance overlap with 35 um. Differences are
   parameter sensitivity, not numerical error.

Formal outputs are written to
`/data/dzf/SecAct_2026/results/spatialess_publication_validation/`.

## Results

All 19 independent processes completed with zero unconfirmed discoveries.

### Reproducibility and resources

Five identical-seed runs produced exactly identical active records,
probabilities, p-values and 0.05 significance calls. Median method time was
0.360 s (CV 1.78%); median complete-process wall time was 3.27 s (CV 4.91%);
median peak RSS was 0.690 GiB (CV 0.28%).

### Permutation resolution

Observed probabilities were exactly identical for 20, 50, 100 and 200
permutations. Relative to the 200-permutation reference, p-value Spearman
correlation increased from 0.966 to 0.989 and 0.995 for 20, 50 and 100
permutations. The corresponding significance-set Jaccard values were 0.940,
0.953 and 0.993, with 9, 7 and 1 threshold disagreements. Method time increased
from 0.350 s at 20 permutations to 2.399 s at 200, while peak RSS remained
0.687-0.695 GiB.

Five 100-permutation seeds all yielded the same 148 significant records and
identical significance sets. P-value Spearman correlation against the first
seed was at least 0.999864. These data support 100 permutations for the
illustrative biological analysis; 20 permutations remain a scaling setting,
not a final inferential recommendation.

### Radius sensitivity

Increasing diffusion radius from 20 to 100 um increased diffusion graph edges
from 28,180 to 641,964. Active-record Jaccard relative to 35 um was
0.979-1.000, and the significance sets were identical in this fixture. For
records active under both radii, probabilities and p-values were identical.
This follows from the integrated binary group-support formulation: radius
changes whether a group pair has spatial support but does not continuously
reweight an already supported group pair. It must not be generalized as a
claim that spatial radius is biologically irrelevant.
