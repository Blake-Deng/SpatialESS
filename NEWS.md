# SpatialESSSPARKLE 0.1.4.9003 (pre-release)

- Bundle SpatialESS 0.1.4.9003 and require that engine version or newer.
- Default inference through the existing H5AD bridge to reference Matrix/BLAS
  group aggregation, retaining strict-tail count/nboot p-values and the explicit
  stream option. Forward aggregation and workspace limits through the existing API.
- Match reference CSR distance rounding and apply ratio after raw-coordinate
  distance calculation; retain the explicit reference and stream choices.
- Include the engine's reference-arithmetic fixes, runtime future.apply dependency,
  regression tests and matched-environment validation summary.
- Automatically check and publish engine/bridge source packages on sparkle-v*
  tags. Development versions are marked as pre-releases; no new biological
  validation or DOI is claimed.

# SpatialESSSPARKLE 0.1.4.9002 (development)

- Bundle the matching SpatialESS Visium HD engine and its CLI/config/validation
  documentation. The existing corrected-H5AD bridge remains a separate entry point.
- Preserve complete branch contents including H5AD manuals and CI.
- Correct installation/issue URLs to the SPARKLE branch of Blake-Deng/SpatialESS.
