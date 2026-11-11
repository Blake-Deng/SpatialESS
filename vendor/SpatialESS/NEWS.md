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
