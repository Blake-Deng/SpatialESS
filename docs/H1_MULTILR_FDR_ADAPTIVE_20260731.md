# Multi-LR FDR and adaptive permutation validation

q^Date: 2026-07-31

> Historical calibration record. The direct/two-sided analytic hybrid described
> below failed the high-permutation gate and is superseded by
> `H1_CONSERVATIVE_HYBRID_VALIDATION_20260731.md`.
^
## Multi-LR BH/FDR calibration

The benchmark contains 10 LR records on a 400-cell, two-group, 16-block
spatial grid:

- 8 simple single-gene LR records;
- 1 two-subunit ligand complex control;
- 1 receptor-cofactor control;
- 2 implanted A-to-B signals in the alternative condition.

The simple records are eligible for the stratified analytic candidate. The
complex and cofactor records are required to use permutation fallback. Each
primary condition uses 499 block-preserving permutations, plus-one correction and the
`greater_equal` tail across 100 null/alternative replicate pairs. A separate
20-replicate reference uses 9,999 permutations per condition.

Results:

| Metric | Analytic simple | Hybrid | Permutation reference |
|---|---:|---:|---:|
| Analytic completion | 80% | 100% after fallback | 100% |
| Null BH FWER | 0.060 | 0.050 | 0.030 |
| Implanted signal power | 1.00 | 1.00 | 1.00 |
| Alternative mean FDP | not applicable to simple-only q | 0.0233 | 0.0233 |
| Simple-LR p Spearman | 0.9664 | - | - |

The 30-replicate pilot is retained in the adjacent result directory, but the
100-replicate run is the primary FDR record. In the 20 x 9,999 reference,
analytic completion remained 80%, implanted-signal power was 1.00 for both
methods, alternative threshold agreement was 99.5%, and mean FDP was 0.0333
for hybrid versus 0.0167 for permutation. Null hybrid FWER was 0.10 versus
0.00 for permutation (2 versus 0 replicates). With only 20 null replicates this
is imprecise, but it prevents enabling analytic/hybrid inference by default.

Machine-readable results:

`benchmarks/results/h1_multilr_fdr_rep100_perm499/`

`benchmarks/results/h1_multilr_fdr_rep20_perm9999/`

## Adaptive permutation screening

`experimental_adaptive_permutation_cellchat_group_support()` uses batches and
simultaneous Clopper-Pearson intervals to stop records whose single-threshold
decision is resolved. The benchmark compares fixed 999 permutations against
`min=100`, `batch=100`, `max=999` on 10 spatial fixtures.

The adaptive run used a mean of 189.95 permutations per record versus 999 for
the fixed reference. Median per-fixture runtime was 0.048 s versus 0.222 s,
and all 10 fixtures had identical 0.05 threshold decisions. Mean absolute p
difference was 0.0220. Because records stop at different sample sizes, these
adaptive p-values are not used for BH/FDR or pathway/network claims.

After the first batch establishes the active domain, the record-level C++
scorer evaluates only records that have not stopped. This used 69.8% of the
full-matrix record-permutation work in the benchmark. Exact sliced-versus-full
probability, reject-count and p-value regression tests use identical seeds.

Machine-readable results:

`benchmarks/results/h1_adaptive_reference/`
