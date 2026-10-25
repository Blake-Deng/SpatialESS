# Reproduce the included results

Install the package, then run `Rscript examples/lr_panel_example.R` for a
database-free example. `Rscript -e 'testthat::test_local(".")'` runs the tests.

The released inference records are in `results/records/`, with LR and
parameter contracts in `results/contracts/`. `results/tables/inputs.tsv`
identifies each prepared input by filename and SHA256. Obtain the external
prepared inputs separately; raw expression matrices are not distributed here.
Their statistical settings are listed in `results/tables/parameters.tsv`.

```sh
Rscript benchmarks/run_prepared.R /path/to/input.rds /path/to/output
```

The runner uses the installed SpatialESS 0.1.4 package, preserves the input's
seed and saves records, parameters, graph sizes, elapsed inference time and
session information. For whole-process peak RSS on Linux:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
/usr/bin/time -v -o run.time.txt \
Rscript benchmarks/run_prepared.R input.rds run_output
```

Single-fixture worker RSS includes input loading and output writing. The
CosMx published-harness resource table uses the complete parent fixture and
is a different measurement scope. Its runner is
`benchmarks/run_secact_cosmx_v3_headtohead.R`; set `SECACT_V3_OUT`,
`SECACT_V3_METHOD=spatialess`, `SECACT_V3_CELLS`, `SECACT_V3_PANEL=full`,
`SECACT_V3_NPERM=20`, and `SECACT_V3_SEED=20260803`. The output directory must
contain the prepared `secact_cosmx_v3_common_fixture.rds`.

Official reference results carry their comparator version and source identity.
Official runtime and RSS measurements come from the completed reference
experiment; they are not a contemporaneous rerun of the October SpatialESS
measurements. Do not present their ratios as a newly paired timing experiment.

The cohort LMM output is the completed 2026-09-07 fit. It uses the current
model implementation and the 45-slice score matrix; fitting date and input
scope are listed separately from the single-sample inference measurements.
It is not labeled as an October LMM refit.

Regenerate the included resource and sparsity figures:

```sh
python -m pip install matplotlib numpy
python scripts/plot_results.py .
```

Verify every distributed file with `sha256sum -c SHA256SUMS.txt` from the
repository root. Environment details are in `results/environment/`.
