# SpatialESS demo

Run a six-cell synthetic example with two measured genes (`L` and `R`), two
cell groups and a seeded permutation analysis. No spatial transcriptomics
input files or GPU are required.

## Install and run

Requires R ≥ 4.1 and a C++17 compiler. Install the package from GitHub:

```r
install.packages(c("remotes", "Matrix", "Rcpp"))
remotes::install_github(
  "Blake-Deng/SpatialESS",
  dependencies = FALSE,
  upgrade = "never"
)
```

From this directory:

```sh
Rscript demo.R
```

The [expected output](expected_output.txt) reports ligand–receptor coverage,
identifies the unmeasured receptor-complex subunit and filters `L_RC`. That
warning is expected in the two-gene example. The retained `L_R` interaction
produces four group-pair records: A→A, A→B, B→A and B→B.

## Measured requirements

Both environments used SpatialESS 0.1.4 with required R dependencies already
installed. Package installation measurements exclude installing R and its
dependencies; the server measurement includes downloading the source.

| Measurement | Apple M3 Pro | Ubuntu server (x86_64) |
| --- | --- | --- |
| R | 4.5.2 | 4.5.3 |
| OS | macOS | Ubuntu 24.04.4 LTS |
| Package installation | Approximately 18 s | 47.70 s |
| Installed package size | Approximately 0.66 MiB | Not recorded |
| Demo whole-process runtime | Median 0.66 s; range 0.50–0.95 s across five runs | Not recorded here |
| Peak resident memory | Approximately 247 MiB | Not recorded here |

The server used Matrix 1.7-5, Rcpp 1.1.1-1.1 and conda-forge GCC 15.2.0.
These measurements describe this small demonstration and do not estimate
performance for the full manuscript analyses.

## Use your own data

Follow the [analysis guide](../docs/QUICKSTART.md) and
[input requirements](../docs/INPUTS.md) for aligned expression, coordinates,
cell groups and ligand–receptor coverage. The
[materials index](../docs/SUBMISSION_CONTENTS.md) links the research analyses.

`demo.R` matches [the package example](../examples/lr_panel_example.R) and is
distributed under [GPL-3](LICENSE).
