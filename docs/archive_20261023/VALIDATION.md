# Candidate validation 2026 10 23

Actual candidate archives were extracted to a newly created system temporary
directory. No old server directory was present. Each figure was run in its own
software extraction with only the unique material objects from its exact input
contract available. Full logs and structured results are delivered separately.

| Check | Execution | Result |
| --- | --- | --- |
| Software and companion SHA-256 | Actually run | Passed |
| Current Figure 1 redraw and composition | Actually run with Figure 1 inputs only | Passed |
| Current Figure 2 redraw | Actually run with Figure 2 inputs only | Passed |
| S1, S2, S3, S4 entries | Actually run independently | All passed |
| Original A fixed-hash negative check | Actual tampered local copy | Rejected before composition |
| R CMD build and install | Actually run from software extraction | Passed |
| Existing testthat suite | Actually run with installed 0.1.4 | Passed |
| Vanilla and instrumented install | Actually run into worker library paths | Both passed |
| All 50,040 feature identities and saved support counts | Actually read saved RDS/TSV and counted nonzero support | Matched Data 4; all-section patient counts include single-section patients |
| TSTA aggregate and cohort sensitivity totals | Actually read saved result tables | Matched manuscript values |
| Core R/C++/Rd/tests/CI | Git content comparison | No changes |
| Original hardware/RCTD historical environment | Static evidence only | Not reconstructed; limitations retained |
| Raw dataset preparation and full LR inference | Not run | Public raw inputs/appropriate compute needed |
| New model fitting, BH or balanced-subset inference | Not run | Saved results retained |
| New Zenodo publication | Not performed | User confirmation still required |

Python environment used: 3.12.14, numpy 2.5.4, pandas 2.2.3, matplotlib 3.11.2,
Pillow 13.0.0, PyMuPDF 1.26.7, openpyxl 3.1.5. R 4.5.2 and a local C++ toolchain
were used. Exact run environments are recorded in new output logs, not presented
as historical environments. Poppler `pdftocairo` was used for Figure 1 composition.
Matplotlib can change PDF metadata and font serialization. The candidate retains
the original final figure objects; newly redrawn figures were validation outputs.

Two known material limitations remain: the copied TSTA_0117 record stream is
truncated; the full-feature balanced-effect summary is absent. Neither is replaced
with invented data. The complete sensitivity score matrix is available as RDS and
can be serialized with `tools/export_saved_scores.R` without inference/model fitting.
