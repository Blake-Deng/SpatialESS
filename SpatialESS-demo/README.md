# SpatialESS demo

This small demo installs SpatialESS 0.1.4 and runs the upstream synthetic ligand–receptor example. It uses six cells, two measured genes (`L` and `R`), two cell groups and one seeded permutation analysis. No spatial transcriptomics input files are needed.

## Tested environment

The demo was tested on the following server environment:

- OS: Ubuntu 24.04.4 LTS (x86_64)
- R: 4.5.3
- SpatialESS: 0.1.4, installed from [`materials-2026.10.23-r1`](https://github.com/Blake-Deng/SpatialESS/tree/materials-2026.10.23-r1)
- Direct R dependencies used by SpatialESS: Matrix 1.7-5 and Rcpp 1.1.1-1.1
- C++ compiler used for the installation: conda-forge GCC 15.2.0 with C++17 support
- GPU: not required; the package and this example run on CPU

The test ran on a server with NVIDIA GPUs available, but neither installation nor the demo used a GPU. Other R versions supported by the package may work, but were not tested for this demo. Installation and runtime below are measurements from the server listed above, not general performance guarantees.

## Desktop validation

The demo was independently tested on an Apple M3 Pro running R 4.5.2. SpatialESS was installed into a temporary local R library and removed after validation; no package installation or demo process remained after testing.

With required R dependencies already available, the temporary SpatialESS installation took approximately 18 s. The installed package occupied approximately 0.66 MiB.

The six-cell synthetic demo was run five times successfully. Whole-process runtime, including R startup, ranged from 0.50 to 0.95 s, with a median of 0.66 s. Peak resident memory was approximately 247 MiB. No GPU or other specialized hardware was required.

These measurements describe this small demonstration only and should not be interpreted as performance estimates for the full spatial transcriptomics analyses reported in the manuscript.

## Install
Measured desktop installation time was approximately 18 s on an Apple M3 Pro with R 4.5.2, with required R dependencies already available. The temporary installed package occupied approximately 0.66 MiB.
Install `remotes` if needed, then install the fixed SpatialESS materials tag:

```r
install.packages(c("remotes", "Matrix", "Rcpp"))
remotes::install_github(
  "Blake-Deng/SpatialESS@materials-2026.10.23-r1",
  dependencies = FALSE,
  upgrade = "never"
)
```

Measured package installation time: **47.70 seconds**. This includes downloading the SpatialESS source from GitHub and compiling/installing the package; it excludes downloading R and R package dependencies. The tested server already had the required dependencies installed. On a new R installation, install the listed dependencies first.

## Run

From this folder, run:

```sh
Rscript demo.R
```

The actual output from a successful run is saved in [`expected_output.txt`](expected_output.txt). It reports ligand–receptor coverage, identifies the missing `unmeasured` receptor-complex subunit and filters `L_RC`. The warning about filtering `L_RC` is expected: only `L` and `R` are present in this two-gene example. The remaining `L_R` interaction produces four group-pair records: A→A, A→B, B→A and B→B. The permutation seed and inference parameters are unchanged from the upstream example.

Measured desktop demo runtime: median 0.66 s across five successful runs (range, 0.50–0.95 s), including R startup, on an Apple M3 Pro with R 4.5.2. Peak resident memory was approximately 247 MiB.

## Use your own data

Replace the example `expression`, `coordinates` and `group` objects with your aligned inputs, then prepare the ligand–receptor table with `filter_cellchat_lr()` and pass the filtered table to `spatialess()`. The input contract covers matrix orientation, identifiers, coordinate alignment, groups, missing genes and database fields: see the main repository's [`docs/INPUTS.md`](https://github.com/Blake-Deng/SpatialESS/blob/materials-2026.10.23-r1/docs/INPUTS.md).

Complete paper code, results and supporting materials are archived at [doi:10.5281/zenodo.23591306](https://doi.org/10.5281/zenodo.23591306).

## License and source

This demo is distributed under GPL-3. `demo.R` is the upstream [`examples/lr_panel_example.R`](https://github.com/Blake-Deng/SpatialESS/blob/materials-2026.10.23-r1/examples/lr_panel_example.R), used without changes. See [`LICENSE`](LICENSE) for the license text and attribution.
