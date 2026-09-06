# SecAct HCC CosMx validation protocol

## Source

- Article: *Inference of secreted protein signaling activities in
  intercellular communication*, DOI `10.1038/s41592-026-03172-0`.
- Data source: Bruker Spatial Biology/NanoString, *CosMx SMI Human Liver
  FFPE Dataset*.
- Official combined Seurat object:
  `LiverDataReleaseSeurat_newUMAP.RDS` (19.7 GiB).
- Official reproduction code: `data2intelligence/SecAct_main`.
- The article's preprocessing code expects a pre-split `CancerousLiver`
  object. `benchmarks/prepare_secact_cosmx_from_bruker.R` reproducibly makes
  this split from the currently distributed combined Seurat object using its
  `Run_Tissue_name` metadata, then applies the published cell-type mapping.

Place the official RDS at
`/data/dzf/SecAct_2026/LiverDataReleaseSeurat_newUMAP.RDS`, then run:

```bash
/home/dzf/miniforge3/envs/cellchat-acceleration/bin/Rscript \
  benchmarks/prepare_secact_cosmx_from_bruker.R
```

The prepared input is written to
`/data/dzf/SecAct_2026/LIHC_CosMx_data.rda`, with source checksum, sample
counts and mapped cell-type counts in adjacent provenance files.

## Fixed primary configuration

- Sample: `CancerousLiver`.
- Expression QC: at least 50 detected genes per cell, matching the SecAct
  reproduction script.
- Normalization: `log1p(count / cell library size * 1000)`, matching the scale
  factor used in the SecAct application.
- Biological groups: the published `metaData$cellType`, not artificial spatial
  grid IDs.
- Spatial graph: binary undirected CSR radius graph with a 20-micrometer
  radius, matching `SecAct.CCC.scST(..., radius = 20)`.
- Molecular score: CellChat-compatible triMean, complex, coreceptor, agonist,
  antagonist and Hill formulas, evaluated only on spatially supported group
  pairs.
- Null: fixed graph, positions and biological groups; complete expression
  profiles are permuted within `sample x niche x 100-micrometer block`.
- Significance: 100 permutations and finite-sample corrected upper-tail
  p-values `(b + 1) / (B + 1)`.

The 100-micrometer block is deliberately larger than the 20-micrometer
communication neighborhood. It retains local tissue and niche structure while
allowing profiles from multiple cell types to be exchangeable. The output
reports the exchangeable-cell fraction so this assumption remains auditable.

## Run

```bash
bash benchmarks/run_secact_cosmx_spatialess.sh
```

For a five-permutation integration check before the primary run:

```bash
SPATIALESS_NPERM=5 \
SPATIALESS_OUT=/data/dzf/SecAct_2026/results/secact_cosmx_smoke \
bash benchmarks/run_secact_cosmx_spatialess.sh
```

`time_verbose.txt` records process-level peak RSS. `summary.tsv` records stage
times, dimensions, graph size, LR coverage and significance counts.

## Interpretation boundary

SecAct estimates downstream secreted-protein activity, whereas SpatialESS
scores ligand-receptor communication on spatial support. Their values are not
expected to be numerically identical. The biological comparison should test
overlap and rank concordance of implicated secreted proteins and
sender-receiver cell-type pairs. Numerical equivalence remains reserved for
the CellChat-compatible reference path.
