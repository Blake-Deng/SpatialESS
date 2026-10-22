# Archived RCTD environment evidence

The archived `spacexr_DESCRIPTION` in this repository identifies version
**2.2.1**, built with R 4.5.3 on 2026-09-03. The saved RCTD script and completed
RAW run log support the use of doublet-mode annotation. The original RCTD run
did not retain `sessionInfo()` or a lockfile, so the complete annotation
environment cannot be reconstructed from the archive.

The currently supplied manuscript wording `spacexr 2.2.16` conflicts with the
archived DESCRIPTION. Use 2.2.1 only if the archived DESCRIPTION is confirmed to
be the package actually loaded by the original annotation process; otherwise
resolve the discrepancy from a contemporaneous package library or run record
before naming an exact version. Do not infer an exact patch version from a
later verification session.

Suggested manuscript wording while the discrepancy is unresolved:

> RCTD doublet-mode labels were taken from an archived annotation run using the
> 17k Ovarian Cancer scFFPE reference and FLEX barcode annotations. The archived
> package record identifies spacexr 2.2.1; the original run did not retain an R
> session record or lockfile, so its complete software environment is
> unavailable.

The saved labels were transferred by `cell_id` and held fixed across raw and
SPARKLE-corrected inputs. This provenance statement does not imply that RCTD
was rerun during the SpatialESS comparison.

