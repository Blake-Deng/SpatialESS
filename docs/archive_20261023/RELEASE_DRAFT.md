# SpatialESS portable figures and reproduction materials 2026 10 23

Proposed tag: `materials-2026.10.23`. Software version: 0.1.4.
Historical release `zenodo-2026.10.22` and DOI 10.5281/zenodo.23567248 are retained.
This candidate is not yet released and has no newly published DOI.

Figure inputs are loaded by target; current Figure 1 composition uses a fixed
published-asset hash instead of server-only paths. Shared results are stored once
in a companion archive with original-to-canonical path/hash mapping. Saved TSTA
feasibility evidence and cohort sensitivity tables/scripts/logs are included.
The R/C++ inference engine, statistical definitions, saved results and conclusions
are unchanged. Inference and model fitting were not rerun.

Actual validation: R build/install and existing regression suite; vanilla and
instrumented source installation; current main figures and S1-S4 from isolated
input sets; fixed Figure 1A tamper rejection; all-feature saved support/identity
checks. Source Data 1/2 and Supplementary Data 1-4 mapping is in the guide.

Known limits: a copied TSTA_0117 records stream is truncated (retained and marked);
the full-feature balanced-effect summary was unavailable (selected summaries,
selection lists and prepared matrix retained). Original RCTD sessionInfo and the
2.2.16-versus-2.2.1 version discrepancy remain unresolved.

Assets to verify before publication: `SpatialESS_software_20261023.tar.gz`,
`SpatialESS_reproduction_20261023.tar.gz`, `FILES.tsv`, `SHA256SUMS.txt`.
If Zenodo's automatic source archive supplies the same software tree, retain one
software copy and manually add/verify the companion. New-version DOI and actual
file inventory must be confirmed before updating citation metadata or manuscript.
