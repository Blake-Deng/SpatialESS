# SPARKLE attribution and scope

SPARKLE is an independent method developed by the SPARKLE authors and
distributed under the MIT License:

- Repository: https://github.com/WangShuai-3/SPARKLE
- Python package: `stambient`
- Validated commit: `e0855e11ce22d4e97fedeeaf5bdcfe94f79b3d78`
- Local version used in the ovarian case study: `0.1.2`
- Preprint: *SPARKLE: evidence-constrained correction of local RNA leakage in
  high-resolution spatial transcriptomics*
- DOI: https://doi.org/10.64898/2026.08.12.744394
- Status at validation: preprint; users should verify the current citation.

This repository does not claim the SPARKLE correction model as a SpatialESS
contribution. It provides a documented bridge that sends corrected expression
to SpatialESS without changing the SpatialCellChat V3-compatible probability
definition. Users must cite SPARKLE when enabling the correction step.

The ovarian benchmark tests computational compatibility and numerical
fidelity after correction. It is not an independent revalidation of every
biological claim made by SPARKLE.
