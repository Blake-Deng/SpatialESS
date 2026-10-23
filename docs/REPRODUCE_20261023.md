# SpatialESS candidate distribution 2026 10 23

This branch prepares a data and plotting revision of SpatialESS 0.1.4. It does
not rerun LR inference, fit new association models, or change scientific results.
The engine, statistical definitions, CI and regression tests are retained.
The historical release `zenodo-2026.10.22` and DOI
https://doi.org/10.5281/zenodo.23567248 remain unchanged. The new companion has
not been published and has no new version DOI. Do not cite the old DOI as the
archive of this candidate's newly added files.

## Downloads and assembly

The two actual candidate files are `SpatialESS_software_20261023.tar.gz` and
`SpatialESS_reproduction_20261023.tar.gz`. The software archive contains the
repository source, real fixture, frozen source snapshots, scripts, configurations,
available environment records and historical execution logs. The companion
contains deduplicated result tables, serialized results, Source Data, figures
and data provenance. Its `objects/` stores one authoritative file per SHA-256;
`MATERIALS.tsv` maps every logical path to that file. The archives do not repeat
this collection. Source aliases in the software tree are relative symlinks to
byte-identical root source or fixtures; the six modified instrumented files
remain independent. Use tar with symlink support on Linux/macOS.

Extract the two archives beside each other. From the software directory:

```sh
python3 tools/verify_software.py
python3 tools/reproduction.py verify --materials ../SpatialESS_reproduction_20261023
python3 -m pip install -r requirements-reproduce.txt
python3 tools/redraw.py Figure1 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py Figure2 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S1 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S2 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S3 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S4 --materials ../SpatialESS_reproduction_20261023
```

`--materials` verifies/mounts only the input set for that figure. The exact sets
are in `archive_20261023/FIGURE_INPUTS.json`. Figure 1 does not load lung models;
Figure 2 does not load CosMx performance or the superseded patient examples.
S4D uses I03/I21/I19; current Figure 2D uses I15/I06/I19. Plots use saved values;
no model/FDR recomputation is part of these commands. `pdftocairo` (Poppler) is
required for Figure 1 vector composition. Its fixed source manifest verifies
the unmodified original A asset and all Figure 1 inputs before/after composition.
The original A SHA-256 is
`92c19a82c3dd3667720db50139d0d043b7d95cad3795d410c54eb5ab5722ec52`.
Arial may be supplied with `SPATIALESS_FONT_DIR`; Liberation Sans is the fallback.
Run exports go to each figure's output directories and `local_runs/`; original
historical output objects are not overwritten. PDF metadata/font rendering may
vary. Frozen inputs and the original A are checked by exact SHA-256.

Mount all material paths for inspection or archived analysis scripts:

```sh
python3 tools/reproduction.py mount --materials ../SpatialESS_reproduction_20261023
```

These links provide compatibility for all historical reading paths, without
copying the data. `MATERIALS.tsv` gives the canonical current storage address,
original size/hash and decoded gzip hash. Equivalent gzip streams are merged
only after comparing unsorted decompressed bytes; barcode order is preserved.
Legacy nested manifests and experiment configs describe the original files and their original compressed hashes. Current verification uses `MATERIALS.tsv`, `SUBMISSION_SHA256.tsv` and the fixed Figure 1 input manifest. `scripts/run_all.py` now dispatches the current portable entries. The old absolute paths inside logs/configuration remain experimental evidence.
Archived preprocessing/inference scripts require public raw inputs and adaptation
of the documented configuration; they are not the portable redraw commands.

## Installation and ablation sources

R >=4.1 and a C++17 toolchain, Matrix and Rcpp are required. R's testthat is
needed for tests; lme4/data.table/future.apply/RANN support optional analyses.

```sh
R CMD build --no-build-vignettes --no-manual .
R CMD INSTALL SpatialESS_0.1.4.tar.gz
Rscript -e 'library(SpatialESS); testthat::test_dir("tests/testthat")'
Rscript tools/install_ablation.R .
```

The ablation installer uses the packaged vanilla and instrumented sources and
installs the exact `lib_vanilla`/`lib` locations consumed by `worker.R`. The
historical design records its original Rscript/server paths. New Linux inference
experiments must select their Rscript and raw fixtures explicitly; the archived
`run_queue.py` uses GNU time, Linux affinity and resource limits. It is not run
for this material revision. `INSTRUMENTED_SOURCE_AUDIT.tsv` and `SOURCE_ALIASES.tsv`
retain implementation identity. `tests/testthat/fixtures/v3_percent_boundary.rds`
is a real regression input and is retained in the software archive.

## Manuscript supplement mapping

| Manuscript item | Actual logical file(s) after mounting | Scope |
| --- | --- | --- |
| Supplementary Data 1 | `paper/submission_materials/Supplementary_Data_1_Ablation.xlsx` | 170 ablation runs and exact summaries |
| Supplementary Data 2 | `paper/submission_materials/Supplementary_Data_2_TSTA.xlsx` | All 1,258 samples, input manifest, status and performance |
| Supplementary Data 3 | `paper/submission_materials/Supplementary_Data_3_Cohort_Sensitivity.xlsx` | Summaries, 28 features, support, influence and all-feature agreement; full model tables are in `paper/cohort_sensitivity/data/` |
| Supplementary Data 4 | `paper/submission_materials/Supplementary_Data_4_Lung_Support.tsv.gz` | All 50,040 primary feature records and zero-inclusive counts |
| Source Data Figure 1 | `paper/submission_materials/Source_Data_Figure1.xlsx` | Supplied final panel workbook; dated figure folder retains its exact workbook too |
| Source Data Figure 2 | `paper/submission_materials/Source_Data_Figure2.zip` | Supplied panel CSVs, identities and counts; workbooks are not assumed interchangeable |
| Supplementary Software 1 | `SpatialESS_software_20261023.tar.gz` | Source, scripts, configuration, database snapshots, contracts, environment evidence and logs; result inputs are supplied in the companion |
| Current figures | `paper/figure_package/updates/Figure1_Atera_20261016/figures/`, `.../Figure2_20261012/figures/` | Original PDF/SVG/PNG exports retained |
| Supplementary figures | `paper/figure_package/figures/FigureS1` through `FigureS4` | Original inputs and entries retained; sensitivity S5 retained separately |

Material files are actual archived objects listed in `MATERIALS.tsv`, not empty
directory promises. `BASELINE_INVENTORY.tsv`, `MIGRATION.tsv` and the added-source
index record original path, size, hash, purpose, known script references and
final location. Filename searches are a conservative reference index, not proof
that an archived script is portable or was executed in this validation.

## TSTA and sensitivity evidence

TSTA sources are at https://doi.org/10.5281/zenodo.17664023. The candidate includes
all-sample summaries, metadata, input manifests, validation reports, scripts,
source and two representative output folders, not all 1,258 result RDS objects.
The aggregate checks give 1,258 completed samples, 132 studies, 37 cancer labels
and 3,648,329 retained spots. This remains SpatialESS-only feasibility evidence;
it is not paired reference validation. TSTA_0117's copied `records.tsv.gz` is
truncated and is retained as received, explicitly excluded as a usable records
table. TSTA_0001's representative records are usable. The historical aggregate
completion record does not certify integrity of every copied representative.

Sensitivity evidence includes complete primary/run/TMA/HC3/Welch tables, fit
states, saved 45-by-50,040 scores (`prepared_slice_scores.rds`), leave-one-patient
summaries, 1,000 balanced patient selections, selected effects, S5 and logs.
The full-feature `balanced_run_adjusted_effects.tsv.gz` and plain score TSV were
not found in the supplied local collection. The matrix is present in the RDS;
no synthetic replacement output or new inferential result was generated.
The historical summary/refit scripts name these additional inputs and are
retained as reconstruction evidence, not claimed to run unchanged end to end.

Descriptive support uses nonzero scores; patient support means at least one
nonzero section. The all-sections count includes single-section patients.
Association tests and support counts remain separate. Zero is not biological
absence; no condition-specific or consensus classification is introduced.
CPU snapshots collected after experiments are not per-run hardware evidence.
The archived spacexr DESCRIPTION says 2.2.1; the manuscript's 2.2.16 and original
RCTD sessionInfo remain unresolved and were not silently changed.

## Publication gate

Proposed release: `materials-2026.10.23` (not created). Both candidate package
hashes and exact sizes are in the delivery `SHA256SUMS.txt`/`FILES.tsv`.
A future fixed URL would be the repository's release download path for this
exact tag and filename; it is not a live download until published.
Use New version of Zenodo record 23567248/concept 23567247. Explicitly upload
and verify the companion: GitHub release attachments are not assumed to be
collected by the webhook. If Zenodo's automatic source snapshot supplies the
software tree adequately, do not upload an identical extra source archive.
After public file-by-file size/hash validation, update the manuscript availability
paragraphs, README, CITATION.cff and .zenodo.json together using the real new
version DOI. This candidate retains the historical DOI metadata pending approval.

The saved section matrix can be exported without inference or fitting:

```sh
Rscript tools/export_saved_scores.R paper/cohort_sensitivity/data/prepared_slice_scores.rds /tmp/slice_scores.tsv.gz
python3 tools/download_public_zenodo.py --record 17664023
```

The second command lists the actual public source files and provider checksums.
Use explicit `--file NAME --out DIRECTORY` to download selected large inputs;
provider checksums are verified and SHA-256 is recorded. Original benchmark input
contracts and manifests remain the authority for the prepared matrices; the
public-source downloader does not establish equivalence of a new preparation.
