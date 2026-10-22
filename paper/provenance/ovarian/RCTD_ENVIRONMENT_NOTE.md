# Archived RCTD environment evidence

The archived `spacexr_DESCRIPTION` in this repository identifies version
**2.2.1**, built with R 4.5.3 on 2026-09-03. The saved RCTD script and completed
RAW run log support the use of doublet-mode annotation. The original RCTD run
did not retain `sessionInfo()` or a lockfile, so the complete annotation
environment cannot be reconstructed from the archive.

An unverified reference elsewhere lists `spacexr 2.2.16`, which conflicts with
the archived DESCRIPTION. Use 2.2.1 only if the archived DESCRIPTION is
confirmed to be the package loaded by the original annotation process;
otherwise resolve the discrepancy from a contemporaneous package library or
run record before naming an exact version. Do not infer an exact patch version
from a later verification session.

The saved labels were transferred by `cell_id` and held fixed across raw and
SPARKLE-corrected inputs. This provenance record does not imply that RCTD was
rerun during the SpatialESS comparison.
