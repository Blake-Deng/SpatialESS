# Figure captions

## Figure 1. Numerical fidelity, performance and spatial-neighborhood sparsity

**A**, Editable summary of sparse inference and patient analysis. The original
workflow asset is preserved separately for provenance. Label permutations reuse
the same spatial graph; group outputs populate the slice-by-feature score matrix.
**B**, Probability agreement for the largest completed paired CosMx fixture
(50,000 cells). Both axes show log10 of positive saved probabilities, with identical
limits and a dashed identity line; hexagon color represents record count. Pearson
correlation is calculated on the original probability scale and maximum absolute
difference retains probability units. Top and right marginal bars count reference
and ESS probabilities in 24 equal log10 bins, retaining every paired record.
**C**, Saved probability correlation, maximum
absolute difference, exact permutation-p fraction, active-record Jaccard and
significance-call disagreement across six completed paired CosMx scales. Tile
colors distinguish metric rows and do not define a shared continuous score;
values rounded to 1 are not a claim of exact floating-point probability equality.
**D,E**, Completed wall-clock runtime and peak RSS measurements. Blue circles
denote SpatialESS and gray squares the official SpatialCellChat V3 comparator.
Logarithmic axes show saved measurements; ESS repeated measurements use their
saved medians. Reference values are archived measurements; the 50,000-cell
annotations report archived reference / current ESS ratios, not newly matched
timing. Vertical guides identify the displayed comparison.
Light fill between completed paired curves is a comparison guide, not a
confidence interval or an extrapolated measurement.
**F-H**, Potential ordered-pair space N squared, observed directed non-self radius
edges E, mean outgoing degree E/N and graph density E/N squared across tested
CosMx, cervical Xenium and ovarian Visium HD fixtures. Only the nested CosMx
ladder is connected in F and H. Panel G shows all 14 fixtures as horizontal
value bars with endpoint symbols and printed degree; its labels identify the
platform, sampling design and cell count. These descriptive bars have no CI.
The light fill in F distinguishes the potential-pair versus observed-edge gap,
not uncertainty. All tested graphs use effective radius 35 micrometres.
Dispersed cluster-stratified Xenium fixtures have lower local degree than central
fixtures; this comparison is descriptive and geometry-dependent.
**I**, Completion and comparison scope. Filled marks indicate paired completion,
open ESS circles indicate ESS-only completion, and crosses indicate recorded
comparator failures. Failures are not treated as completed numerical performance
or fidelity observations. CosMx significance disagreement uses the saved raw
permutation-p threshold of 0.05, not a BH-adjusted discovery claim.

## Figure 2. Patient-aware condition associations in lung fibrosis

**A**, The frozen primary patient-random-intercept LMM analysis, with PF minus
Control effect on saved log1p communication scores and minus log10 of saved
BH q. All valid fits are shown; convergence warnings, numerical failures and
insufficient-nonzero fits remain outside the primary discovery set. The original
eight labeled features show full LR names and retain their stable IDs. Every
label has its own vector leader terminating at the saved full-key effect and
BH-q coordinate, with an arrowhead and a highlighted target. **B**, The
original 12 selected effects as zero-to-estimate horizontal bars, overlaid with
saved 95% Wald CIs and estimate markers. The aligned right-hand column prints
the same effects in units of 0.001; open points denote valid singular fits and
filled points non-singular fits. **C**, The original 28 communication features, with
ligand/receptor, stable feature ID and sender/receiver identity. All displayed
rows have the same sender and receiver group, which does not imply same-cell
self-communication. LMM and slice-OLS effect / response-SD values share one
zero-centered color scale and are printed in their tiles. Dots indicate each
model's saved feature-level BH q below 0.05. Support is the nonzero-slice fraction
with denominator 45; patient ICC and singular status have separate encodings.
Complex receptor displays follow the audited molecular annotation, including
CSPG4-(ITGA3+ITGB1), without changing the original feature key or statistics.
**D**, The three original examples, identified by LR, stable ID and actual
sender-to-receiver groups. Small points show all 45 slices per feature; outlined
points show each of the 35 patient means of saved log1p scores. Fixed horizontal
jitter changes only display positions. Condition-effect CIs for these examples
appear in Supplementary Figure S4D. No adjusted group-mean CI is inferred.
**E**, The original 52 pathway/cell-pair representatives. Color shows saved
standardized effect, and area is proportional to minus log10 of the saved LR
feature's BH q, capped at 10. These values are not pathway-enrichment FDR.
Blank positions have no saved representative. Complete feature identities,
unabbreviated groups and values are in the source-data workbook.

## Supplementary Figure S1. Existing SPARKLE bridge validation

**A**, The documented corrected-H5AD input bridge and comparison on identical
inputs. **B**, Twelve archived engine-regression or bridge/direct comparisons from
six input configurations (RAW/corrected at three sizes);
E and B identify these distinct comparisons. **C**, Archived official/SpatialESS
corrected-input probability agreement with density bins, an identity line and
reference/ESS marginal record-count distributions in 24 equal log10 bins.
**D**, Permutation-p agreement; point area is proportional to the number of
records sharing each p-value pair. **E**, Active, exact-p and raw-p interaction-call
agreement, displayed as horizontal value bars with endpoint markers and values.
**F**, Shared, lost and gained raw-p significant records when comparing
RAW and SPARKLE-corrected inputs. Panel F measures sensitivity to changed inputs,
not numerical engine equivalence. Its significance threshold is raw p below
0.05. Archived runs are explicitly identified and are not presented as newly
rerun official measurements.

## Supplementary Figure S2. Independent cervical Xenium validation

**A**, Actual coordinates for the 5,000-cell cluster-stratified fixture, with the
original five most frequent public graph-cluster labels and remaining clusters
grouped as Other. Equal spatial aspect and coordinate units are preserved.
**B**, Current saved fidelity metrics for eight central/stratified size and filter
configurations, with the documented dimensions-only comparator repair. Row
suffixes indicate min.percent / min.cells.sr. **C**, Full-key paired probabilities
for four default-filter fixtures; both axes use log10, matched limits and an
identity line, with marginal count distributions from the same paired records.
**D**, Saved active and raw-p significant record counts as paired horizontal
bars and numeric columns by fixture and filter, blue circles for ESS and gray
squares for V3. The threshold is raw permutation p
below 0.05; no BH-FDR interpretation or biological cell-type annotation is added.

## Supplementary Figure S3. Array storage and sampling geometry

**A**, Theoretical dense float64 N-by-N storage (8N squared bytes; open points)
and measured CSR offsets, neighbors and distances (filled points) for the five
original representative fixtures. The logarithmic axis uses decimal MB. These
quantities are not measured whole-process peak RSS. **B**, Mean local degree in
central and cluster-stratified cervical fixtures, shown as horizontal value bars.
**C**, Fraction of cells with zero outgoing radius-graph degree in all 14 tested
fixtures, displayed as horizontal bars with printed percentages and platform
symbols. Tiny nonzero percentages retain sufficient precision in their labels.
**D**, CosMx measured
core-CSR array payload and all stored numeric-array payload, including the
additional graph fields listed in source data; light fill marks the measured
gap between these array payloads, not a confidence interval. Radius, input and directed-edge
definitions match Figure 1F-H; no graph or communication inference is rerun.

## Supplementary Figure S4. Model diagnostics and full group labels

**A**, Mutually exclusive saved LMM fit states as horizontal bars, markers and
aligned counts. Valid singular fits are distinguished
from convergence warnings and numerical failures. **B**, Saved BH-q discovery
overlap between the patient LMM and slice-OLS analysis, using the same value-bar
style. **C**, The same 52 bubble
representatives as Figure 2E, with expanded sender/receiver mappings. Color and
size definitions and limits are identical to Figure 2E. **D**, Saved condition
effects as zero-to-estimate bars and saved 95% Wald CIs for the three Figure 2D
examples. These intervals describe
the PF minus Control condition contrast, not adjusted condition-group means.
