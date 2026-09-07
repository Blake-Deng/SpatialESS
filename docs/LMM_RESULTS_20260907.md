# GSE250346 robust LMM: completed results

This is a results update for the distributed 2026-09-07 robust-final code.
No core model source was changed while producing these tables. The cohort
contains 45 slices from 35 patients: 9 Control and 26 PF patients.
See [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md) for interpretation and limitations.

## Result index

Compact TSVs are in [results/tables](../results/tables/):

| File | Content |
|---|---|
| `gse250346_lmm_model_summary.tsv` | Cohort, formula, fit counts and parameter lock |
| `gse250346_lmm_fit_status.tsv` | Complete mutually exclusive status counts |
| `gse250346_glm_vs_lmm_summary.tsv` | Significant counts, overlap and effect agreement |
| `gse250346_lmm_before_after_summary.tsv` | Baseline versus robust comparison |
| `gse250346_lmm_status_transitions.tsv` | Exact status transitions |
| `gse250346_lmm_convergence_summary.tsv` | Warning classes with documented precedence |
| `gse250346_lmm_failure_summary.tsv` | Errors remaining after all three optimizers |
| `gse250346_lmm_optimizer_summary.tsv` | Optimizer sequences attempted |
| `gse250346_lmm_non_singular_sensitivity.tsv` | Nonsingular-only analysis; both BH universes labeled |
| `gse250346_lmm_support_diagnostics.tsv` | Median support measures by status |
| `gse250346_lmm_result_audit.tsv` | Eight checks, all passed |

[gse250346_lmm_robust_20260907](../results/tables/gse250346_lmm_robust_20260907/)
contains full LMM and GLM feature tables (`.tsv.gz`), per-feature transitions,
sample metadata, input checksums, source snapshot, session information and
run/finalization logs. [gse250346_lmm_baseline](../results/tables/gse250346_lmm_baseline/)
preserves the previous summaries and normalized full feature tables.
No expression matrices or patient-identifying information are bundled.

All 50,040 features are present. Valid fits number 39,570; 1,984 have warnings,
292 have numerical failures and 8,194 have insufficient nonzero support.
BH-adjusted q-values are NA for every non-`ok` feature. There are 5,147
significant LMM features using the stated Wald-normal p-values and q < 0.05.
These are model-based screening associations, not proven disease mechanisms.

## Frozen computation

```text
LMM: log1p(score) ~ condition + (1 | patient)
GLM: log1p(score) ~ condition
Contrast: conditionDisease (PF versus Control)
REML: FALSE
min_units: 6
min_nonzero: 2
max_zero_fraction: 0.98
min_response_sd: 1e-12
rank_tolerance: 1e-10
optimizers: bobyqa, nloptwrap, Nelder_Mead
p_method: Wald_normal
BH: all valid features together, not separately within worker chunks
workers: 32 feature processes; 1 native/BLAS thread per process
R: 4.5.3; lme4: 2.0-1; Matrix: 1.7-5
```

The source file `R/multisample.R` is identical in SpatialESS, the SPARKLE
vendor copy and the previously supplied robust-final code packages:

```text
SHA256: 8cd1241b2a052d07a8882cf04787bd8bcf846966057dab5a226fe86911cf20fa
MD5:    5bffbc01532d1b8c1f0ee49cbfbcf449
```

The final table was merged from 1,001 completed feature checkpoints. Its
`elapsed_seconds` is deliberately NA on resume: `invocation_elapsed_seconds`
measures finalization only, not model fitting. These tables are not a new
SpatialESS-versus-official runtime or peak-memory benchmark. A metadata-only
change to the number of workers does not mean communication inference was
rerun. This update fits downstream cohort models from the existing 45 inputs.

## Read and verify without rerunning fits

From the repository root, using only base R:

```bash
Rscript benchmarks/verify_gse250346_result_bundle.R results/tables
```

```r
con <- gzfile("results/tables/gse250346_lmm_robust_20260907/patient_random_intercept_lmm.tsv.gz", "rt")
tab <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
close(con)
stopifnot(nrow(tab) == 50040L)
table(tab$status)
hits <- subset(tab, status == "ok" & !is.na(q_value) & q_value < 0.05)
```

TSVs use quoting because optimizer messages may contain embedded newlines.
Use a structured TSV reader, not line counting or splitting on newlines.
Baseline tables were normalized from the original unquoted multiline format
by the supplied legacy reader; numerical fields and original statuses are
preserved. Original source files remain untouched on the benchmark server.

## Reproduce fits from prepared inputs

The input root requires `full_sample_results/*.rds` (45 SpatialESS output
profiles), `sample_metadata_from_rds.tsv` and the original baseline folder
`slice_glm_vs_patient_lmm/` for historical comparison. These large prepared
inputs are not in this table update; hashes are supplied for provenance.
See [REPRODUCE_MAIN.md](REPRODUCE_MAIN.md) for the upstream workflow.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript benchmarks/run_gse250346_lmm_checkpointed.R \
  . /path/to/multisample /path/to/new_output 4
Rscript benchmarks/audit_gse250346_lmm_tables.R \
  /path/to/new_output /path/to/multisample/slice_glm_vs_patient_lmm
```

Use up to the available CPU allocation, not more workers than available cores.
The completed run used 32 workers. A checkpoint signature checks source,
inputs, R and lme4 before resuming; use a fresh directory if they change.

GEO accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346
This cohort did not undergo SPARKLE correction. Its inclusion in the
SPARKLE repository validates the common downstream engine only.
