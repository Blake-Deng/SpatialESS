test_that("v3-compatible streaming matches a direct edge reference", {
  expression <- Matrix::Matrix(
    matrix(
      c(1, 0.5, 0,
        0, 1, 0.5),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("L", "R"), paste0("c", 1:3))
    ),
    sparse = TRUE
  )
  coordinates <- matrix(
    c(0, 0, 1, 0, 3, 0), ncol = 2, byrow = TRUE,
    dimnames = list(colnames(expression), c("x", "y"))
  )
  group <- factor(c("A", "A", "B"))
  lr <- data.frame(
    ligand = "L", receptor = "R", interaction_name = "L_R",
    pathway_name = "test", annotation = "Secreted Signaling",
    agonist = "", antagonist = "", co_A_receptor = "", co_I_receptor = "",
    row.names = "L_R", stringsAsFactors = FALSE
  )
  complex <- data.frame(subunit_1 = character(), row.names = character())
  cofactor <- data.frame(cofactor1 = character(), row.names = character())

  result <- spatial_v3_compatible(
    expression, coordinates, group, lr, complex, cofactor,
    ratio = 1, tol = 0, interaction.range = 3,
    scale.distance = 1, contact.range = 1,
    min.percent = 0, min.cells.sr = 1,
    nboot = 0
  )

  ligand <- as.numeric(expression["L", ])
  receptor <- as.numeric(expression["R", ])
  distance <- as.matrix(stats::dist(coordinates))
  minimum <- min(distance[distance > 0])
  probability <- matrix(0, 3, 3)
  for (sender in seq_len(3)) {
    for (receiver in seq_len(3)) {
      if (sender != receiver && distance[sender, receiver] > 3) next
      product <- ligand[sender] * receptor[receiver]
      if (product <= 0) next
      spatial <- if (sender == receiver) 1 / minimum else 1 / distance[sender, receiver]
      probability[sender, receiver] <- product / (0.5 + product) * spatial
    }
  }
  expected <- list()
  index <- 1L
  for (sender_group in levels(group)) {
    for (receiver_group in levels(group)) {
      values <- probability[group == sender_group, group == receiver_group, drop = FALSE]
      values <- values[values > 0]
      if (!length(values)) next
      expected[[index]] <- data.frame(
        sender_group_name = sender_group,
        receiver_group_name = receiver_group,
        probability = mean(values)
      )
      index <- index + 1L
    }
  }
  expected <- do.call(rbind, expected)
  observed <- result$records[, c(
    "sender_group_name", "receiver_group_name", "probability"
  )]
  key <- function(x) paste(x$sender_group_name, x$receiver_group_name)
  observed <- observed[match(key(expected), key(observed)), , drop = FALSE]
  expect_identical(
    observed[, c("sender_group_name", "receiver_group_name")],
    expected[, c("sender_group_name", "receiver_group_name")]
  )
  expect_equal(observed$probability, expected$probability, tolerance = 1e-15)
})
