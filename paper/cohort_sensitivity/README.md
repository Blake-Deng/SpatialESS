# Cohort model sensitivity

The saved analysis includes 50,040-feature model tables, fit states, the
45-section score RDS, patient-omission summaries and execution logs. Download
and mount the companion data package using the
[reproduction guide](../../docs/REPRODUCE.md).

The 4,904/5,021 retained primary discoveries and 3,808 descriptive cross-analysis
support count match the saved sensitivity tables. Selected balanced summaries
and all 1,000 patient selections are retained. The all-feature balanced-effect
table is unavailable. Scripts include their input contracts and recorded paths;
full reruns require the stated inputs and environment configuration.

The [condition-support analysis](../lung_cohort/condition_support/README.md)
provides a score TSV, descriptive support labels and one-section-per-patient
checks. These labels are distinct from the model sensitivity results and do
not introduce new condition-association tests or BH correction.
