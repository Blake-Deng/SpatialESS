# ROI numerical-audit scripts

These are the scripts used for the reported ROI preparation, original-RCTD
comparison, and result QC. They do not use SSH or download data.

1. Create a fresh RUN_DIR with `inputs/official_outs` pointing to the verified
   official Space Ranger output and `inputs/ovarian_reference.rds` containing
   `list(counts=gene_by_cell_raw_UMI, cell_types=named_labels)`.
2. Run `Rscript spatialess_prepare_real_hd_roi.R RUN_DIR`. The default shared ROI
   is 5×5 complete 32µm blocks with >=80% tissue/observed/depth coverage. The reported
   ROI started at 32µm grid row 100, column 104 and had 100% coverage.
3. The generated audit configs use 5×5 blocks. For the reported bin8 numerical
   audit only, change its `roi32$n_rows` and `$n_cols` to 2 before running it;
   leave the production CLI config at 5×5. Run
   `Rscript spatialess_real_rctd_audit.R CONFIG_RDS` for each 8/16/32 config.
4. Run the repository CLI with `inputs/real_roi_cli_config.rds` and then
   `Rscript spatialess_real_roi_qc.R RUN_DIR`.

Read each script's `--help`. Paths are supplied as arguments/config, and outputs
are protected against unintended overwrite. All settings and selected geometry
are recorded. This audit checks numerical fidelity, not independent biology.
These scripts target the original spacexr 2.2.1 API and the stated SpatialESS
development version; they are not a promise of compatibility with other APIs.
