make_integrated_fixture <- function() {
  cell_id <- paste0("c", seq_len(80))
  position <- rep(seq_len(20), 4)
  block <- rep(seq_len(4), each = 20)
  group <- factor(ifelse(position %% 2L, "A", "B"), levels = c("A", "B"))
  expression <- Matrix::Matrix(
    rbind(
      Ld = 1 + (position %% 5), Rd = 1 + ((position + 2) %% 5),
      Lc = 1 + (position %% 3), Rc = 1 + ((position + 1) %% 3)
    ),
    sparse = TRUE,
    dimnames = list(c("Ld", "Rd", "Lc", "Rc"), cell_id)
  )
  coordinates <- cbind(x = position, y = block * 100)
  rownames(coordinates) <- cell_id
  lr <- data.frame(
    interaction_name = c("Ld_Rd", "Lc_Rc"),
    pathway_name = c("diff", "contact"),
    ligand = c("Ld", "Lc"), receptor = c("Rd", "Rc"),
    annotation = c("Secreted Signaling", "Cell-Cell Contact"),
    stringsAsFactors = FALSE
  )
  complex <- data.frame(subunit_1 = character(), row.names = character())
  cofactor <- data.frame(cofactor1 = character(), row.names = character())
  list(expression = expression, coordinates = coordinates, group = group,
       lr = lr, complex = complex, cofactor = cofactor)
}

test_that("integrated workflow uses mechanism-specific CSR graphs", {
  fixture <- make_integrated_fixture()
  result <- spatialess(
    fixture$expression, fixture$coordinates, fixture$group,
    fixture$lr, fixture$complex, fixture$cofactor,
    contact_radius = 2, diffusion_radius = 8, block_size = 100,
    inference = "permutation", nperm = 9, seed = 71,
    verbose = FALSE
  )
  expect_s3_class(result, "SpatialESSIntegrated")
  expect_setequal(names(result$graphs), c("contact", "diffusion"))
  expect_equal(result$graphs$contact$parameters$radius, 2)
  expect_equal(result$graphs$diffusion$parameters$radius, 8)
  expect_lt(length(result$graphs$contact$neighbors),
            length(result$graphs$diffusion$neighbors))
  expect_setequal(result$records$mechanism, c("contact", "diffusion"))
  expect_true(all(result$records$permutation_confirmed))
  expect_true(all(result$records$inference_mode ==
                  "fixed_spatial_permutation"))
  expect_equal(result$records$interaction_name,
               result$lr$interaction_name[result$records$global_lr_index])
})

test_that("guarded integrated discoveries are permutation-confirmed", {
  fixture <- make_integrated_fixture()
  result <- spatialess(
    fixture$expression, fixture$coordinates, fixture$group,
    fixture$lr, fixture$complex, fixture$cofactor,
    contact_radius = 2, diffusion_radius = 8, block_size = 100,
    inference = "guarded_hybrid", nperm = 19, seed = 83,
    alpha = 0.05, screening_margin = 0.01, verbose = FALSE
  )
  significant <- result$records$pvalue <= result$parameters$alpha
  expect_true(all(result$records$permutation_confirmed[significant]))
  expect_identical(result$parameters$inference, "guarded_hybrid")
  expect_true(all(vapply(
    result$inference,
    function(x) identical(x$parameters$mode, "guarded_hybrid"),
    logical(1)
  )))
})

test_that("integrated workflow enforces aligned cells", {
  fixture <- make_integrated_fixture()
  bad_coordinates <- fixture$coordinates[
    rev(seq_len(nrow(fixture$coordinates))), , drop = FALSE
  ]
  expect_error(
    spatialess(
      fixture$expression, bad_coordinates, fixture$group,
      fixture$lr, fixture$complex, fixture$cofactor,
      nperm = 2, verbose = FALSE
    ),
    "cell order differs"
  )
})
