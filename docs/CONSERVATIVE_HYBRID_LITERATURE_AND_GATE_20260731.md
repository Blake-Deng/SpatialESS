# Conservative analytic screening with permutation confirmation

Date: 2026-07-31

## Decision

Stage 1B uses the analytic approximation only to exclude clearly
non-significant records.

For target level `alpha` and non-negative margin `m`:

1. If a record is eligible for the guarded analytic model and
   `p_analytic > alpha + m`, set its hybrid p-value to 1 and classify it as
   analytically screened non-significant.
2. Otherwise, compute the fixed block-preserving permutation p-value.
3. Complex, cofactor, agonist, antagonist, autocrine, zero-quartile,
   heavy-tie, small-group, insufficient-shared-block and high-dependence
   records always use permutation.
4. A hybrid result at or below `alpha` must have
   `permutation_confirmed = TRUE`.
5. BH/FDR and publication p-values use fixed permutation. The adaptive engine
   remains a single-threshold experimental accelerator and is not used to
   report BH-adjusted p-values.

The implementation is
`experimental_hybrid_cellchat_group_support()` in `R/hybrid.R`.

## Why this is conservative

Let `P_i` be the p-value from the complete fixed-permutation reference for
record `i`. Define the hybrid value as

`P_i^H = P_i` for records sent to permutation, and `P_i^H = 1` for
analytically screened records.

Therefore `P_i^H >= P_i` for every record and every realized data set.
Any coordinatewise monotone multiple-testing rule, including Bonferroni,
Holm, BH and BY, can only produce a subset of the discoveries produced from
the full permutation vector. The screen cannot introduce a discovery or
increase empirical FWER/FDR relative to the same full-permutation reference.

This argument does not prove that the full-permutation BH procedure controls
FDR under every pattern of dependent LR tests. BH has its own independence or
positive-dependence assumptions; BY is available for arbitrary dependence,
and resampling-based maxT/step-down procedures require their own conditions.
The benchmark must consequently report both absolute null error and direct
hybrid-versus-permutation agreement.

The cost of the one-way rule is possible false-negative screening. Exact
agreement is not automatic. Stage 1B requires the hybrid and full-permutation
rejection sets to be identical over the calibration matrix before the screen
can be enabled beyond experimental use.

## Fixed versus adaptive permutation

Random-permutation publication p-values use the plus-one estimator

`(b + 1) / (B + 1)`,

where `b` is the number of permuted statistics at least as extreme as the
observed statistic and `B` is the number of random permutations. This avoids
zero Monte Carlo p-values and is the convention used in the formal fixed
reference.

Sequential Monte Carlo can save work near obvious decisions, but optional
stopping requires explicit resampling-risk control. The current adaptive
prototype uses simultaneous confidence bounds over planned looks for a
single threshold. It is not treated as an exact p-value generator and is not
used for BH/FDR. Anytime-valid Monte Carlo p-values or an MMCTest-style
multiple-testing procedure are candidates for a later, separately validated
adaptive implementation.

## Literature search

The search used Crossref and OpenAlex metadata APIs and targeted queries for
sequential Monte Carlo tests, adaptive Monte Carlo multiple testing, random
permutation p-values, independent filtering, and anytime-valid Monte Carlo
p-values. The following papers directly inform the design.

1. Besag J, Clifford P (1991). Sequential Monte Carlo p-values.
   Biometrika. https://doi.org/10.1093/biomet/78.2.301
   Establishes sequential Monte Carlo testing based on ordered stopping
   rules. It supports adaptive computation, not unrestricted reuse of a
   naively stopped exceedance fraction as a fixed-sample p-value.

2. Gandy A (2009). Sequential implementation of Monte Carlo tests with
   uniformly bounded resampling risk. Journal of the American Statistical
   Association. https://doi.org/10.1198/jasa.2009.tm08368
   Formalizes the probability that a Monte Carlo implementation makes a
   different decision from the ideal test. This is the relevant safety target
   for adaptive threshold decisions.

3. Sandve GK, Ferkingstad E, Nygard S (2011). Sequential Monte Carlo multiple
   testing. Bioinformatics.
   https://doi.org/10.1093/bioinformatics/btr568
   Extends adaptive allocation to multiple hypotheses and makes clear that
   multiplicity and simulation uncertainty must be handled jointly.

4. Gandy A, Hahn G (2014). MMCTest: a safe algorithm for implementing
   multiple Monte Carlo tests. Scandinavian Journal of Statistics.
   https://doi.org/10.1111/sjos.12085
   Provides a multiple-testing framework with an explicit probability of
   disagreement with the ideal multiple-testing result.

5. Gandy A, Hahn G (2017). QuickMMCTest: quick multiple Monte Carlo testing.
   Statistics and Computing.
   https://doi.org/10.1007/s11222-016-9656-z
   Prioritizes simulation effort across uncertain hypotheses. It is relevant
   for future batching, but does not justify the current fixed-reference FDR
   results by itself.

6. Fischer L, Ramdas A (2024 preprint; 2026 journal record). Multiple testing
   with anytime-valid Monte Carlo p-values.
   https://doi.org/10.48550/arxiv.2404.15586
   Develops p-values valid under continuous monitoring and optional stopping.
   This is a stronger future basis for adaptive BH than a naively stopped
   Monte Carlo fraction.

7. Phipson B, Smyth GK (2010). Permutation p-values should never be zero:
   calculating exact p-values when permutations are randomly drawn.
   https://doi.org/10.2202/1544-6115.1585
   Supports the plus-one correction used by the fixed reference.

8. Hemerik J, Goeman J (2018). Exact testing with random permutations.
   TEST. https://doi.org/10.1007/s11749-017-0571-1
   Clarifies when random-permutation tests retain finite-sample validity and
   the role of the identity transformation and exchangeability.

9. Bourgon R, Gentleman R, Huber W (2010). Independent filtering increases
   detection power for high-throughput experiments. PNAS.
   https://doi.org/10.1073/pnas.0914005107
   Independent filtering can preserve type-I error when the filter is
   independent of the test statistic under the null. SpatialESS does not
   assume that condition for analytic p-values, so this theorem is not used
   to justify the screen.

10. Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate.
    Journal of the Royal Statistical Society B.
    https://doi.org/10.1111/j.2517-6161.1995.tb02031

11. Benjamini Y, Yekutieli D (2001). The control of the false discovery rate
    in multiple testing under dependency. Annals of Statistics.
    https://doi.org/10.1214/aos/1013699998
    These two papers delimit what can be claimed from BH under dependent LR
    tests and motivate a BY sensitivity analysis.

12. Westfall PH, Young SS (1993). Resampling-Based Multiple Testing.
    Wiley. https://doi.org/10.1002/9781118057516
    Provides the reference framework for resampling-based family-wise error
    control; strong-control claims require assumptions such as subset
    pivotality.

## Stage 1B acceptance criteria

- Run 100 null/alternative replicates with 499 fixed permutations.
- Run 20 null/alternative replicates with 9,999 fixed permutations.
- Expand group size, zero fraction, ties and LR dependence.
- Assert every hybrid discovery is permutation-confirmed.
- Assert hybrid p-values are coordinatewise no smaller than the matching full
  permutation values.
- Require exact hybrid/full-permutation rejection-set agreement throughout
  the declared calibration domain.
- Require null FWER/FDR to be no greater than the matching permutation
  reference; report confidence intervals rather than interpreting a small
  number of repetitions as proof.
- Keep hybrid experimental and do not freeze Stage 1B if any exact-agreement
  or error-control gate fails.
