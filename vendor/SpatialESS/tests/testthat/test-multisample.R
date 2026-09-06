make_sample_records <- function(scores) {
  feature <- c("F1", "F2", "F3")
  data.frame(
    sender_group_name = "Sender",
    receiver_group_name = "Receiver",
    interaction_name = feature,
    probability = unname(scores),
    stringsAsFactors = FALSE
  )
}

test_that("patient-level preparation aggregates slices before modelling", {
  sample_id <- paste0("s", seq_len(8))
  patient_id <- rep(paste0("p", seq_len(4)), each = 2)
  condition <- rep(c("control", "control", "disease", "disease"), each = 2)
  sample_metadata <- data.frame(
    sample_id = sample_id,
    patient_id = patient_id,
    condition = condition,
    stringsAsFactors = FALSE
  )
  scores <- list(
    c(1, 1, 5), c(1, 1, 5),
    c(1, 1, 5), c(1, 1, 5),
    c(1, 5, 1), c(1, 5, 1),
    c(1, 5, 1), c(1, 5, 1)
  )
  names(scores) <- sample_id
  results <- lapply(scores, make_sample_records)
  names(results) <- sample_id

  prepared <- prepare_multisample_communication(
    results, sample_metadata,
    unit = "patient"
  )
  expect_s3_class(prepared, "SpatialESSMultiSample")
  expect_equal(nrow(prepared$score), 4)
  f2 <- prepared$feature_metadata$feature_id[
    prepared$feature_metadata$interaction_name == "F2"
  ]
  expect_equal(prepared$score["p1", f2], 1)
  expect_equal(prepared$score["p3", f2], 5)

  model <- fit_multisample_communication_glm(
    prepared,
    design = ~ condition,
    coefficient = "conditiondisease",
    min_units = 4,
    min_nonzero = 1
  )
  f2 <- model[model$interaction_name == "F2", , drop = FALSE]
  f3 <- model[model$interaction_name == "F3", , drop = FALSE]
  expect_equal(f2$status, "ok")
  expect_equal(f3$status, "ok")
  expect_gt(f2$estimate, 0)
  expect_lt(f3$estimate, 0)

  conserved <- summarize_conserved_communication(
    prepared,
    min_prevalence = 1,
    min_units_per_group = 2
  )
  expect_true(conserved$conserved[conserved$interaction_name == "F1"])
})

test_that("patient metadata must be constant within patient", {
  sample_metadata <- data.frame(
    sample_id = c("s1", "s2"),
    patient_id = c("p1", "p1"),
    condition = c("control", "disease"),
    stringsAsFactors = FALSE
  )
  results <- list(
    s1 = make_sample_records(c(1, 1, 1)),
    s2 = make_sample_records(c(1, 1, 1))
  )
  expect_error(
    prepare_multisample_communication(results, sample_metadata),
    "varies within patient"
  )
})
