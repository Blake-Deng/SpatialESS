# SpatialESS algorithm specification v0

Date: 2026-07-30

Status: design baseline; formulas and decision gates must be validated before
the public API is frozen.

## 1. Scientific target

SpatialESS will infer ligand-receptor communication only over biologically
plausible spatial support and will quantify uncertainty under a null model that
preserves tissue structure. Its practical cost should scale with sparse spatial
edges and supported LR edges, rather than an explicit all-cell `N x N` domain.
No universal linear-complexity claim will be made because output size, graph
density, LR count and solver iterations remain data-dependent.

## 2. Inputs and invariants

Required inputs:

- non-negative gene-by-cell expression matrix, preferably `dgCMatrix`;
- cell identifiers and 2D/3D coordinates;
- cell state/type labels;
- LR, complex and cofactor definitions;
- optional sample, patient, tissue compartment and segmentation identifiers.

Invariants:

- cell identifiers must align exactly across expression, coordinates and
  metadata;
- graph indices and edge counts use 64-bit-safe representations;
- self-edges, duplicate edges and cross-sample edges are explicit policy
  choices, never silent side effects;
- contact and diffusible interactions use distinct graphs;
- all random procedures expose and record a seed;
- sparse matrices stay sparse unless a bounded conversion is documented.

## 3. Two scoring paths

### 3.1 Compatibility path: SpatialESS-G

This path is the controlled baseline. It retains CellChat type-7 triMean,
complex/cofactor handling and Hill activation. Spatial support determines which
directed sender-receiver group pairs are eligible and supplies a documented
spatial factor. It provides a bridge to CellChat/SpatialCellChat and a numerical
reference for the new method.

For LR pair `l` and sender/receiver groups `(a,b)`, let `L_la` and `R_lb` be the
CellChat-compatible ligand and receptor summaries. The molecular component is

`M_lab = (L_la R_lb)^n / (Kh^n + (L_la R_lb)^n)`.

Only group pairs supported by at least one mechanism-appropriate spatial edge
are materialized. The exact definition of the spatial factor must be selected
by benchmark; candidates are supported-edge fraction, normalized edge mass and
distance-decay mass. Raw molecular and spatial components will both be retained
in output so their effects are auditable.

### 3.2 New mechanistic path: SpatialESS-F

This path is a new model and will not be described as numerically identical to
CellChat. Contact LR pairs are evaluated on segmentation/radius adjacency.
Diffusible ligands generate a concentration field on a sparse tissue graph.

Let `W` be a non-negative spatial weight matrix and `L_G = D - W` its graph
Laplacian. For ligand `l`, the initial linear reaction-diffusion model is

`dc_l/dt = alpha_l s_l - (delta_l I + D_l L_G + U_l diag(r_l)) c_l`,

where `s_l` is ligand production support, `r_l` is receptor availability,
`D_l` is graph diffusion, `delta_l > 0` is decay and `U_l >= 0` is optional
receptor-mediated uptake.

At steady state,

`(delta_l I + D_l L_G + U_l diag(r_l)) c_l = alpha_l s_l`.

For a symmetric non-negative graph this system is sparse positive definite and
can be solved with conjugate gradient plus preconditioning. A solve is required
per ligand, not per cell pair. Receiver activation for receptor `r` is then

`A_lri = (c_li R_ri)^n / (Kh^n + (c_li R_ri)^n)`.

Sender attribution will first be reported at cell-state or niche level using a
small number of right-hand sides. Per-cell source attribution will not be
materialized for million-cell data unless explicitly requested.

The dynamic ODE trajectory is optional. The steady-state sparse solve is the
primary scalable method; time integration is reserved for time-series data,
RNA-velocity-informed experiments or deliberate transient simulations.

## 4. Mechanism-specific graph construction

Four initial graph builders are required:

1. segmentation contact adjacency;
2. fixed-radius graph;
3. adaptive k-nearest-neighbor graph with a maximum-radius guard;
4. distance-decay graph with optional compartment barriers.

Every graph returns a canonical edge table containing `sample_id`, `sender`,
`receiver`, `distance`, `weight`, `graph_type` and barrier/compartment flags,
plus a sparse matrix view. Graph construction must never create cross-sample
edges.

## 5. Edge-first execution

For LR pair `l`, the evaluated domain is

`E_l = E_mechanism intersect ligand-supported senders intersect receptor-supported receivers`.

The implementation will intersect sorted sparse index lists rather than scan
all cell pairs. Primary outputs are chunked LR-annotated edge tables or sparse
tensors. Dense group-by-group arrays are derived outputs for compatibility and
visualization only.

## 6. Statistical calibration

The initial reference null is a graph-preserving, sample-preserving permutation
that shuffles labels or expression assignments within tissue compartments or
spatial blocks. Global cell shuffling is not the default.

Planned inference modes:

- `permutation`: high-permutation reference and fallback;
- `analytic`: only for LR/graph classes with demonstrated type-I calibration;
- `hybrid`: analytic screening followed by adaptive spatial permutation for
  borderline, sparse, dependent or complex/cofactor cases.

FastCCC's distribution convolution is a statistical inspiration, not a score
replacement. CellChat triMean/Hill definitions and graph dependence require
independent derivation and calibration.

## 7. Borrowed concepts and method boundaries

- **CellChatESS:** sparse referenced-gene execution, type-7 triMean,
  complex/cofactor/Hill definitions and deterministic threading.
- **kidney interaction framework:** sparse neighborhood construction and
  compartment-preserving label permutation.
- **CytoSignal:** separate contact/diffusion neighborhoods, cell-resolved
  spatial outputs and spatially local signaling analyses.
- **FastCCC:** distribution-based analytic nulls, convolution/FFT and hybrid
  significance motivation.

No external method's score will be relabeled as SpatialESS. Reused code, if
any, must be license-compatible and explicitly attributed; the preferred route
is an independent implementation from documented formulas.

## 8. Algorithm iterations and decision gates

### G0: graph correctness

- exact agreement with brute-force distances on synthetic 2D/3D coordinates;
- no missing/duplicate/cross-sample edges;
- deterministic output under cell reordering;
- memory proportional to stored edges.

### G1: static edge-first score

- exact molecular component versus CellChatESS on supported group pairs;
- exact agreement with a small brute-force spatial implementation;
- explicit contact versus diffusion behavior;
- measurable savings as supported-edge fraction decreases.

### F1: steady-state graph ODE

- solver residual below a declared tolerance;
- non-negative concentration under valid inputs;
- mass-balance and limiting-case tests;
- `D -> 0` and barrier tests behave as predicted;
- robustness curves for radius, diffusion, decay and uptake parameters.

### H1: calibrated significance

- type-I error and FDR on null simulations;
- power on implanted spatial LR signals;
- agreement with 10,000-100,000 spatial permutations on tractable datasets;
- stratification by group size, zero fraction, complex/cofactor status,
  autocrine/paracrine class and spatial dependence.

### Integrated release

- comparison with CellChat/SpatialCellChat, CytoSignal and other selected
  spatial CCC methods;
- one small biological concordance dataset;
- one semi-simulated ground-truth benchmark;
- one million-cell Xenium/CosMx/MERFISH stress test;
- runtime, peak RSS, completion, calibration and biological recovery reported
  separately.

## 9. Immediate implementation order

1. Implement and test graph object plus radius/kNN/contact constructors.
2. Build a brute-force oracle for small coordinate sets.
3. Implement expression-support intersection and chunked edge output.
4. Add the static compatibility score.
5. Add the sparse steady-state ODE solver and limiting-case tests.
6. Add compartment/block-preserving permutation.
7. Only then prototype analytic/hybrid significance.

