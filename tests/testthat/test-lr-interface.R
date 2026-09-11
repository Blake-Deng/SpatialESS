test_that("H5AD filtering matches direct inference and records missing subunits", {
  python <- Sys.getenv("PYTHON", unset = Sys.which("python3"))
  skip_if(!nzchar(python), "Python is unavailable")
  available <- suppressWarnings(system2(python,
    c("-c", shQuote("import anndata, numpy, pandas, scipy")), stdout = FALSE, stderr = FALSE))
  skip_if(available != 0L, "Python bridge dependencies are unavailable")
  root <- tempfile("bridge_lr_"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  script <- file.path(root, "make.py"); h5ad <- file.path(root, "input.h5ad")
  writeLines(c(
    "import anndata as ad, numpy as np, pandas as pd, sys",
    "from scipy import sparse",
    "x = sparse.csr_matrix(np.array([[3,1,100],[1,2,100],[2,3,100],[4,2,100],[2,4,100],[1,1,100]],dtype=float))",
    "obs = pd.DataFrame({'group':['A','A','A','B','B','B'],'x':np.arange(6),'y':np.zeros(6)},index=['c'+str(i) for i in range(6)])",
    "ad.AnnData(x,obs=obs,var=pd.DataFrame(index=['L','R','EXTRA'])).write_h5ad(sys.argv[1])"
  ), script)
  expect_equal(system2(python, shQuote(c(script, h5ad))), 0L)
  lr <- data.frame(ligand = c("L", "L"), receptor = c("R", "RC"),
    interaction_name = c("L_R", "L_RC"), annotation = "Secreted Signaling")
  complex <- data.frame(subunit_1 = "R", subunit_2 = "MISSING", row.names = "RC")
  cofactor <- data.frame(cofactor1 = character(), row.names = character())
  args <- list(h5ad = h5ad, group_col = "group", coordinate_cols = c("x", "y"),
    lr = lr, complex = complex, cofactor = cofactor, gene_list = c("L", "R", "EXTRA"),
    python = python, output_dir = file.path(root, "bundle"),
    nboot = 19L, seed.use = 7L, scale.distance = 1,
    min.percent = 0, min.cells.sr = 1)
  expect_warning(result <- do.call(spatialess_from_h5ad, c(args, list(lr_missing = "filter"))),
                 "Excluded 1 of 2")
  rng <- .Random.seed
  expect_equal(result$bridge$excluded_lr, 1)
  expect_identical(result$lr_filter$audit$missing_receptor_genes, c("", "MISSING"))
  expect_true(file.exists(file.path(root, "bundle", "lr_missing_genes.tsv")))
  b <- read_spatialess_bundle(file.path(root, "bundle"))
  direct <- SpatialESS::spatialess(b$expression, b$coordinates, b$group,
    lr = lr[1, , drop = FALSE], complex = complex, cofactor = cofactor,
    nboot = 19L, seed.use = 7L, scale.distance = 1, min.percent = 0, min.cells.sr = 1)
  expect_identical(result$records, direct$records)
  expect_identical(rng, .Random.seed)
  expect_error(do.call(spatialess_from_h5ad, c(args, list(lr_missing = "error"))), "MISSING")
  expect_warning(legacy <- do.call(spatialess_from_h5ad, args), "Excluded 1 of 2")
  expect_identical(legacy$records, result$records)
  expect_error(do.call(spatialess_from_h5ad, c(args, list(filter_unresolved = FALSE, lr_missing = "filter"))),
               "conflicting policies")
  broken <- args; broken$complex <- data.frame(wrong_column = "R")
  expect_error(do.call(spatialess_from_h5ad, broken), "subunit")
})
