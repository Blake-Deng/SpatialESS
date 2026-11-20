percent_filter_fixture <- function(group, ligand, receptor) {
  expression <- Matrix::Matrix(rbind(L = ligand, R = receptor), sparse = TRUE)
  colnames(expression) <- paste0("cell", seq_along(group))
  list(expression = expression, coordinates = cbind(seq_along(group), 0),
       group = group, lr = data.frame(ligand = "L", receptor = "R"),
       complex = data.frame(subunit_1 = character()),
       cofactor = data.frame(cofactor1 = character()),
       tol = 0, interaction.range = 1, scale.distance = 1, min.cells.sr = 1L)
}

test_that("the default cutoff excludes rare expression written in scientific notation", {
  old <- options(scipen = 0, OutDec = ".")
  on.exit(options(old), add = TRUE)
  group <- factor(rep(c("rare", "tenth", "all"), c(10000L, 10L, 1L)))
  expression <- as.integer(!duplicated(group))
  args <- percent_filter_fixture(group, expression, expression)
  # Each group has one positive cell. R formats the rare fraction as "1e-04".
  numeric <- do.call(spatialess, args)
  expect_identical(numeric$parameters$min.percent, .1)
  expect_identical(numeric$parameters$percent_filter, "numeric")
  expect_setequal(numeric$records$sender_group_name, c("tenth", "all"))
  expect_setequal(numeric$records$receiver_group_name, c("tenth", "all"))

  legacy <- do.call(spatialess, c(args, list(percent_filter = "spatialcellchat_v3")))
  expect_true(any(legacy$records$sender_group_name == "rare" &
                  legacy$records$receiver_group_name == "rare"))

  options(scipen = 999, OutDec = ",")
  explicit <- do.call(spatialess, c(args, list(percent_filter = "numeric")))
  expect_identical(explicit$records, numeric$records)
})

test_that("observed and permuted groups use the same selected expression cutoff", {
  old <- options(scipen = 0, OutDec = ".")
  on.exit(options(old), add = TRUE)
  group <- factor(c(rep(c("A", "B", "C"), 10), "A", "B"))
  within_group <- ave(seq_along(group), group, FUN = seq_along)
  # Ligand in A and receptor in B are positive in 5/11 cells. Both fractions
  # in C are exactly 5/10. A and B have a higher positive ligand dose than C.
  ligand <- (within_group <= c(A = 5, B = 7, C = 5)[group]) *
    ifelse(group == "C", .5, 1)
  receptor <- as.integer(within_group <= c(A = 7, B = 5, C = 5)[group])
  args <- c(percent_filter_fixture(group, ligand, receptor),
            list(min.percent = .5, nboot = 25L, seed.use = 7L))

  # C -> A has observed score .5. All active edges have unit distance weight,
  # so a permutation exceeds .5 exactly when it passes both group cutoffs and
  # a high-dose C sender touches a positive A receiver.
  set.seed(args$seed.use)
  invisible(future.apply::future_lapply(1L, identity, future.seed = TRUE))
  exceeds <- replicate(args$nboot, {
    labels <- group[sample.int(length(group))]
    senders <- which(labels == "C" & ligand == 1)
    receivers <- which(labels == "A" & receptor > 0)
    higher <- any(abs(outer(senders, receivers, "-")) <= 1)
    fractions <- data.frame(ligand = tapply(ligand > 0, labels, mean),
                            receptor = tapply(receptor > 0, labels, mean))
    legacy_mask <- format(fractions, digits = 1) >= .5
    c(numeric = higher && fractions["C", "ligand"] >= .5 &&
        fractions["A", "receptor"] >= .5,
      spatialcellchat_v3 = higher && legacy_mask["C", "ligand"] &&
        legacy_mask["A", "receptor"])
  })
  expected_p <- rowMeans(exceeds)
  expect_lt(expected_p[["numeric"]], expected_p[["spatialcellchat_v3"]])
  for (aggregation in c("reference", "stream")) {
    for (mode in names(expected_p)) {
      result <- do.call(spatialess, c(args, list(aggregation = aggregation, percent_filter = mode)))
      records <- result$records
      expect_identical(any(records$sender_group_name == "A"), mode == "spatialcellchat_v3")
      expect_identical(any(records$receiver_group_name == "B"), mode == "spatialcellchat_v3")
      ca <- records$sender_group_name == "C" & records$receiver_group_name == "A"
      expect_identical(records$probability[ca], .5)
      expect_identical(records$pvalue[ca], expected_p[[mode]])
    }
  }
})
