.feature_key <- function(sender, receiver, interaction_name) {
  paste(
    nchar(sender), sender,
    nchar(receiver), receiver,
    nchar(interaction_name), interaction_name,
    sep = ":"
  )
}

.result_records <- function(result, sample_id) {
  records <- if (inherits(result, "SpatialESSResult")) {
    result$records
  } else if (is.data.frame(result)) {
    result
  } else {
    stop(
      sprintf(
        "Result '%s' must be a SpatialESSResult or data frame.",
        sample_id
      ),
      call. = FALSE
    )
  }
  required <- c(
    "sender_group_name", "receiver_group_name",
    "interaction_name", "probability"
  )
  if (!all(required %in% colnames(records))) {
    stop(
      sprintf("Result '%s' is missing required record columns.", sample_id),
      call. = FALSE
    )
  }
  records <- as.data.frame(records, stringsAsFactors = FALSE)
  records$sender_group_name <- as.character(records$sender_group_name)
  records$receiver_group_name <- as.character(records$receiver_group_name)
  records$interaction_name <- as.character(records$interaction_name)
  records$probability <- as.double(records$probability)
  if (anyNA(records[, required[seq_len(3L)], drop = FALSE]) ||
      any(!is.finite(records$probability)) ||
      any(records$probability < 0)) {
    stop(sprintf("Result '%s' contains invalid records.", sample_id),
         call. = FALSE)
  }
  records$feature_id <- .feature_key(
    records$sender_group_name,
    records$receiver_group_name,
    records$interaction_name
  )
  if (anyDuplicated(records$feature_id)) {
    stop(
      sprintf("Result '%s' contains duplicate communication features.", sample_id),
      call. = FALSE
    )
  }
  records$sample_id <- sample_id
  records
}

#' Prepare patient- or sample-level communication scores
#'
#' Combines independently inferred SpatialESS results into a common
#' sample-by-communication-feature matrix. Missing features are represented by
#' zero. When patient-level analysis is requested, multiple slices from one
#' patient are aggregated before statistical modelling, so cells are never
#' treated as biological replicates.
#'
#' @param results Named list of SpatialESSResult objects or record data
#'   frames.
#' @param sample_metadata Data frame containing one row per input sample.
#' @param sample_col Column containing names that match names(results).
#' @param patient_col Patient identifier column.
#' @param condition_col Biological condition column.
#' @param covariates Additional patient/sample-level covariates to retain.
#' @param unit Whether the output rows represent patients or samples.
#' @param replicate_aggregation Aggregation used for multiple slices from one
#'   patient.
#' @param fill Score assigned when a communication feature is absent.
#' @return A SpatialESSMultiSample object.
#' @export
prepare_multisample_communication <- function(
    results,
    sample_metadata,
    sample_col = "sample_id",
    patient_col = "patient_id",
    condition_col = "condition",
    covariates = character(),
    unit = c("patient", "sample"),
    replicate_aggregation = c("mean", "median"),
    fill = 0) {
  unit <- match.arg(unit)
  replicate_aggregation <- match.arg(replicate_aggregation)
  if (!is.list(results) || !length(results) ||
      is.null(names(results)) || any(!nzchar(names(results))) ||
      anyDuplicated(names(results))) {
    stop("results must be a non-empty uniquely named list.", call. = FALSE)
  }
  if (!is.data.frame(sample_metadata)) {
    stop("sample_metadata must be a data frame.", call. = FALSE)
  }
  required_metadata <- unique(c(
    sample_col,
    if (unit == "patient") patient_col,
    condition_col,
    covariates
  ))
  if (!all(required_metadata %in% colnames(sample_metadata))) {
    stop("sample_metadata is missing required columns.", call. = FALSE)
  }
  if (anyNA(sample_metadata[, required_metadata, drop = FALSE])) {
    stop("Required sample metadata contains missing values.", call. = FALSE)
  }
  sample_id <- as.character(sample_metadata[[sample_col]])
  if (anyDuplicated(sample_id)) {
    stop("sample_metadata contains duplicate sample identifiers.", call. = FALSE)
  }
  missing_metadata <- setdiff(names(results), sample_id)
  if (length(missing_metadata)) {
    stop(
      sprintf(
        "Missing metadata for result samples: %s",
        paste(missing_metadata, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  sample_metadata <- sample_metadata[
    match(names(results), sample_id), , drop = FALSE
  ]
  rownames(sample_metadata) <- names(results)

  records <- Map(.result_records, results, names(results))
  records <- do.call(rbind, records)
  rownames(records) <- NULL
  feature_id <- unique(records$feature_id)
  sample_count <- length(results)
  score <- matrix(
    as.double(fill),
    nrow = sample_count,
    ncol = length(feature_id),
    dimnames = list(names(results), feature_id)
  )
  row_index <- match(records$sample_id, rownames(score))
  column_index <- match(records$feature_id, colnames(score))
  score[cbind(row_index, column_index)] <- records$probability

  first_record <- match(feature_id, records$feature_id)
  optional <- intersect(
    c("pathway_name", "annotation", "ligand", "receptor"),
    colnames(records)
  )
  feature_metadata <- records[first_record, c(
    "feature_id", "sender_group_name", "receiver_group_name",
    "interaction_name", optional
  ), drop = FALSE]
  rownames(feature_metadata) <- feature_metadata$feature_id

  retained_metadata <- unique(c(
    if (unit == "patient") patient_col else sample_col,
    condition_col,
    covariates
  ))
  if (unit == "sample") {
    unit_metadata <- sample_metadata[, retained_metadata, drop = FALSE]
    rownames(unit_metadata) <- as.character(unit_metadata[[sample_col]])
    output_score <- score[rownames(unit_metadata), , drop = FALSE]
  } else {
    patient <- as.character(sample_metadata[[patient_col]])
    patient_levels <- unique(patient)
    constant_columns <- unique(c(condition_col, covariates))
    for (patient_id in patient_levels) {
      selected <- patient == patient_id
      for (column in constant_columns) {
        values <- unique(as.character(sample_metadata[[column]][selected]))
        if (length(values) != 1L) {
          stop(
            sprintf(
              "Metadata column '%s' varies within patient '%s'.",
              column, patient_id
            ),
            call. = FALSE
          )
        }
      }
    }

    output_score <- matrix(
      0,
      nrow = length(patient_levels),
      ncol = ncol(score),
      dimnames = list(patient_levels, colnames(score))
    )
    for (index in seq_along(patient_levels)) {
      selected <- patient == patient_levels[[index]]
      values <- score[selected, , drop = FALSE]
      output_score[index, ] <- if (replicate_aggregation == "mean") {
        colMeans(values)
      } else {
        apply(values, 2L, stats::median)
      }
    }
    first_sample <- match(patient_levels, patient)
    unit_metadata <- sample_metadata[
      first_sample, retained_metadata, drop = FALSE
    ]
    rownames(unit_metadata) <- patient_levels
    unit_metadata[[patient_col]] <- patient_levels
  }

  structure(
    list(
      score = output_score,
      feature_metadata = feature_metadata,
      unit_metadata = unit_metadata,
      sample_metadata = sample_metadata,
      parameters = list(
        unit = unit,
        sample_col = sample_col,
        patient_col = patient_col,
        condition_col = condition_col,
        covariates = covariates,
        replicate_aggregation = replicate_aggregation,
        fill = fill,
        input_samples = sample_count
      )
    ),
    class = "SpatialESSMultiSample"
  )
}

#' Fit sample-level differential communication models
#'
#' Fits one generalized linear model per sender-receiver-LR feature. The rows
#' of prepared$score, not individual cells, are the statistical replicates.
#'
#' @param prepared A result from prepare_multisample_communication().
#' @param design A right-hand-side model formula, for example
#'   ~ condition + age + sex.
#' @param coefficient Name or index of the model coefficient to test.
#' @param family Response model. gaussian_log1p includes zeros and is the
#'   default for continuous communication probabilities; binomial_presence
#'   tests detection prevalence.
#' @param min_units Minimum complete patient/sample units.
#' @param min_nonzero Minimum units with positive communication score.
#' @param presence_threshold Threshold used by binomial_presence.
#' @param adjust_method Multiple-testing correction method.
#' @return Feature-level coefficient, uncertainty and adjusted p-value table.
#' @export
fit_multisample_communication_glm <- function(
    prepared,
    design,
    coefficient,
    family = c("gaussian_log1p", "gaussian", "binomial_presence"),
    min_units = 6L,
    min_nonzero = 2L,
    presence_threshold = 0,
    adjust_method = "BH") {
  if (!inherits(prepared, "SpatialESSMultiSample")) {
    stop("prepared must be a SpatialESSMultiSample object.", call. = FALSE)
  }
  if (!inherits(design, "formula") || length(design) != 2L) {
    stop("design must be a right-hand-side formula.", call. = FALSE)
  }
  family <- match.arg(family)
  min_units <- as.integer(min_units)
  min_nonzero <- as.integer(min_nonzero)
  if (is.na(min_units) || min_units < 3L ||
      is.na(min_nonzero) || min_nonzero < 0L) {
    stop("min_units or min_nonzero is invalid.", call. = FALSE)
  }

  model_frame <- stats::model.frame(
    design, data = prepared$unit_metadata, na.action = stats::na.fail
  )
  model_matrix <- stats::model.matrix(design, data = model_frame)
  if (nrow(model_matrix) != nrow(prepared$score)) {
    stop("Design rows do not match communication score rows.", call. = FALSE)
  }
  coefficient_index <- if (is.character(coefficient)) {
    match(coefficient, colnames(model_matrix))
  } else {
    as.integer(coefficient)
  }
  if (length(coefficient_index) != 1L || is.na(coefficient_index) ||
      coefficient_index < 1L || coefficient_index > ncol(model_matrix)) {
    stop("coefficient is absent from the design matrix.", call. = FALSE)
  }
  coefficient_name <- colnames(model_matrix)[coefficient_index]

  score <- prepared$score
  result <- data.frame(
    feature_id = colnames(score),
    coefficient = coefficient_name,
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    nonzero_units = colSums(score > presence_threshold),
    total_units = nrow(score),
    status = "not_fitted",
    stringsAsFactors = FALSE
  )
  if (nrow(score) < min_units) {
    result$status <- "insufficient_units"
  } else {
    response_family <- if (family == "binomial_presence") {
      stats::binomial()
    } else {
      stats::gaussian()
    }
    for (feature in seq_len(ncol(score))) {
      values <- score[, feature]
      if (result$nonzero_units[feature] < min_nonzero) {
        result$status[feature] <- "insufficient_nonzero"
        next
      }
      response <- switch(
        family,
        gaussian_log1p = log1p(values),
        gaussian = values,
        binomial_presence = as.double(values > presence_threshold)
      )
      if (!is.finite(stats::var(response)) || stats::var(response) == 0) {
        result$status[feature] <- "constant_response"
        next
      }
      fit <- tryCatch(
        stats::glm.fit(
          x = model_matrix,
          y = response,
          family = response_family
        ),
        error = function(error) error
      )
      if (inherits(fit, "error") || !isTRUE(fit$converged)) {
        result$status[feature] <- "fit_failed"
        next
      }
      class(fit) <- c("glm", "lm")
      coefficients <- summary(fit)$coefficients
      if (coefficient_index > nrow(coefficients) ||
          anyNA(coefficients[coefficient_index, ])) {
        result$status[feature] <- "rank_deficient"
        next
      }
      result$estimate[feature] <- coefficients[coefficient_index, 1L]
      result$std_error[feature] <- coefficients[coefficient_index, 2L]
      result$statistic[feature] <- coefficients[coefficient_index, 3L]
      result$p_value[feature] <- coefficients[coefficient_index, 4L]
      result$status[feature] <- "ok"
    }
  }
  result$q_value <- stats::p.adjust(result$p_value, method = adjust_method)
  result <- merge(
    prepared$feature_metadata,
    result,
    by = "feature_id",
    all.y = TRUE,
    sort = FALSE
  )
  result[match(colnames(score), result$feature_id), , drop = FALSE]
}

#' Summarize communication programs conserved across conditions
#'
#' @param prepared A result from prepare_multisample_communication().
#' @param group_col Metadata column defining biological groups.
#' @param presence_threshold Minimum score counted as present.
#' @param min_prevalence Minimum prevalence required in every group.
#' @param min_units_per_group Minimum patient/sample units per group.
#' @return Feature table with group prevalence, group means and a conserved flag.
#' @export
summarize_conserved_communication <- function(
    prepared,
    group_col = NULL,
    presence_threshold = 0,
    min_prevalence = 0.5,
    min_units_per_group = 2L) {
  if (!inherits(prepared, "SpatialESSMultiSample")) {
    stop("prepared must be a SpatialESSMultiSample object.", call. = FALSE)
  }
  if (is.null(group_col)) group_col <- prepared$parameters$condition_col
  if (!group_col %in% colnames(prepared$unit_metadata)) {
    stop("group_col is absent from unit metadata.", call. = FALSE)
  }
  if (!is.finite(min_prevalence) ||
      min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be between zero and one.", call. = FALSE)
  }

  group <- as.character(prepared$unit_metadata[[group_col]])
  group_levels <- unique(group)
  group_sizes <- table(factor(group, levels = group_levels))
  if (any(group_sizes < as.integer(min_units_per_group))) {
    stop("At least one group has too few patient/sample units.", call. = FALSE)
  }
  score <- prepared$score
  prevalence <- mean_score <- matrix(
    0,
    nrow = ncol(score),
    ncol = length(group_levels),
    dimnames = list(colnames(score), make.names(group_levels, unique = TRUE))
  )
  for (index in seq_along(group_levels)) {
    selected <- group == group_levels[[index]]
    prevalence[, index] <- colMeans(
      score[selected, , drop = FALSE] > presence_threshold
    )
    mean_score[, index] <- colMeans(score[selected, , drop = FALSE])
  }
  colnames(prevalence) <- paste0("prevalence_", colnames(prevalence))
  colnames(mean_score) <- paste0("mean_", colnames(mean_score))
  summary <- cbind(
    prepared$feature_metadata,
    as.data.frame(prevalence),
    as.data.frame(mean_score)
  )
  summary$minimum_group_prevalence <- apply(prevalence, 1L, min)
  summary$conserved <- summary$minimum_group_prevalence >= min_prevalence
  rownames(summary) <- NULL
  summary
}
