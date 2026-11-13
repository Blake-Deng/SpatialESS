# SpatialESS 0.1.4.9003 (pre-release)

- Default to reference Matrix/BLAS group aggregation to preserve strict-tail
  counts at exact floating-point ties. Keep the old stream aggregation only as
  an explicit option, and guard dense cell-by-group/group-by-group workspaces.
- Match single-gene, named-complex, R mean/quantile/prod and arithmetic order
  in the CellChat-compatible scorers; restore unambiguous cached components.
- Preserve separate double products in Euclidean distances and apply coordinate
  ratio after raw-distance calculation, matching the reference order.
- Declare future.apply as a runtime dependency and always reproduce the tested
  V3 RNG advance before generating label permutations.
- Keep strict greater-than comparisons and uncorrected count/nboot p-values.
  No epsilon, p-value replacement or change to the reference equations is used.
- Add exact regression tests and matched-environment server validation: 17,880
  active records across six real and five synthetic cases match fresh official
  probabilities, p-values and permutation sequences exactly.
- Add tag-triggered source-package releases after package checks; four-component
  development versions are automatically marked as pre-releases.

# SpatialESS 0.1.4.9002 (development)

- Add Visium HD native sparse readers, raw-UMI grid aggregation, and independent
  RCTD full/doublet/multi fits for requested8/16/32µm bins.
- Preserve filtered feature-ID policy before symbol mapping, complete tissue
  positions, zero-read positions, per-bin exclusions and upstream classifications.
- Fit global RCTD platform effects/sigma once per resolution before pixel batches.
- Add production CLI, per-resolution config, native-count validator, documentation
  and tests; retain RNA-contribution and dominant-bin communication boundaries.
- Record real-data server engineering validation with explicit ROI/native32 and
  biological-validation limits. No formal release or new DOI is asserted here.
- Add a dedicated cloud Visium HD job with pinned original spacexr and explicit
  H5/Parquet dependencies; require zero skipped HD tests. Keep generic package
  checks separate from optional backend resolution.
