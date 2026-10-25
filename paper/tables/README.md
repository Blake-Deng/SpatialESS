# Supplementary Table 1

The TSV/CSV/XLSX files are supplied through the companion data package.
Follow the [reproduction guide](../../docs/REPRODUCE.md) to download and mount it.

The canonical table has 327 rows and 52 fields. Rows are
measurements/comparisons, not independent biological datasets. `benchmark_type`
and `measurement_scope` explicitly distinguish SpatialESS single-fixture repetitions,
CosMx full-parent harness summaries, official fidelity, failed references, the H5AD
bridge and archived ovarian compatibility. TSV/CSV preserve full numerical precision;
XLSX provides the same table plus definitions, environment and CosMx repetitions.
Failed attempt durations and RSS are recorded only in failure-attempt fields.
Dense theoretical storage is not substituted for peak RSS. NA means unavailable,
unmeasured or not applicable. All input rows are sourced from included tables;
model results retain their saved fitting dates and performance measurements retain
their recorded execution scopes.
