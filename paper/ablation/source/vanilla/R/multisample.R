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
    if (unit == "patient") patient_col else c(sample_col, patient_col),
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

#' Fit patient-random-intercept linear mixed models
#'
#' Fits one linear mixed model per communication feature while retaining one
#' row per spatial slice. The response is log1p(probability), the fixed
#' effects are supplied by `design`, and `patient_col` defines a random
#' intercept for patient-level correlation among slices.
#'
#' The returned p-values are Wald normal-approximation p-values. They are
#' intended for a transparent, lightweight screening analysis; confirmatory
#' analyses may use lmerTest or parametric bootstrap degrees of freedom.
#'
#' @param prepared A result from prepare_multisample_communication(...,
#'   unit = "sample").
#' @param design A right-hand-side fixed-effect formula, for example
#'   `~ condition + tma`.
#' @param coefficient Name or index of the fixed-effect coefficient to test.
#' @param patient_col Patient identifier column in prepared$unit_metadata.
#' @param min_units Minimum number of slice-level units.
#' @param min_nonzero Minimum units with positive communication score.
#' @param max_zero_fraction Maximum allowed fraction of zero scores. Set to 1
#'   to disable this support screen while still recording zero_fraction.
#' @param min_response_sd Minimum standard deviation on the log1p-score scale.
#'   Features below this threshold are marked insufficient_variation.
#' @param rank_tolerance Tolerance passed to the fixed-effect design rank check.
#' @param adjust_method Multiple-testing correction method.
#' @param REML Whether to use restricted maximum likelihood. Fixed-effect
#'   comparisons generally use `FALSE`.
#' @param standardize_response Whether to standardize each feature response
#'   before fitting. Estimates and standard errors are returned on the original
#'   log1p-score scale; this can improve conditioning for very small scores.
#' @param optimizer lme4 optimizer used for each feature.
#' @param optimizer_fallbacks Additional optimizers tried after `optimizer`.
#'   A successful warning-free fit is preferred; otherwise the successful fit
#'   with the first available optimizer is retained and flagged.
#' @return Feature-level mixed-model estimates, standard errors, Wald
#'   statistics and p-values, plus optimizer, variance-component,
#'   convergence, singular-fit and failure diagnostics, including the attempted
#'   optimizer sequence. BH-adjusted p-values are reported only for fits with
#'   status equal to ok.
#' @export
fit_multisample_communication_lmm <- function(
    prepared,
    design,
    coefficient,
    patient_col = NULL,
    min_units = 6L,
    min_nonzero = 2L,
    max_zero_fraction = 0.98,
    min_response_sd = 1e-12,
    rank_tolerance = 1e-10,
    adjust_method = "BH",
    REML = FALSE,
    standardize_response = TRUE,
    optimizer = "bobyqa",
    optimizer_fallbacks = c("nloptwrap", "Nelder_Mead")) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for the mixed model.", call. = FALSE)
  }
  if (!inherits(prepared, "SpatialESSMultiSample") ||
      !identical(prepared$parameters$unit, "sample")) {
    stop("prepared must be a sample-level SpatialESSMultiSample object.",
         call. = FALSE)
  }
  if (!inherits(design, "formula") || length(design) != 2L) {
    stop("design must be a right-hand-side formula.", call. = FALSE)
  }
  if (is.null(patient_col)) patient_col <- prepared$parameters$patient_col
  if (!is.character(patient_col) || length(patient_col) != 1L ||
      !patient_col %in% colnames(prepared$unit_metadata)) {
    stop("patient_col is absent from sample metadata.", call. = FALSE)
  }
  min_units <- as.integer(min_units)
  min_nonzero <- as.integer(min_nonzero)
  if (is.na(min_units) || min_units < 3L ||
      is.na(min_nonzero) || min_nonzero < 0L) {
    stop("min_units or min_nonzero is invalid.", call. = FALSE)
  }
  if (!is.numeric(max_zero_fraction) || length(max_zero_fraction) != 1L ||
      !is.finite(max_zero_fraction) || max_zero_fraction < 0 ||
      max_zero_fraction > 1) {
    stop("max_zero_fraction must be a number between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(min_response_sd) || length(min_response_sd) != 1L ||
      !is.finite(min_response_sd) || min_response_sd < 0) {
    stop("min_response_sd must be a non-negative finite number.",
         call. = FALSE)
  }
  if (!is.numeric(rank_tolerance) || length(rank_tolerance) != 1L ||
      !is.finite(rank_tolerance) || rank_tolerance <= 0) {
    stop("rank_tolerance must be a positive finite number.", call. = FALSE)
  }
  if (!is.logical(standardize_response) || length(standardize_response) != 1L ||
      is.na(standardize_response)) {
    stop("standardize_response must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(optimizer) || length(optimizer) != 1L ||
      !nzchar(optimizer)) {
    stop("optimizer must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.character(optimizer_fallbacks) ||
      any(!nzchar(optimizer_fallbacks))) {
    stop("optimizer_fallbacks must contain non-empty names.", call. = FALSE)
  }

  data <- prepared$unit_metadata
  data$.patient_random_effect <- factor(data[[patient_col]])
  model_frame <- stats::model.frame(
    design, data = data, na.action = stats::na.fail
  )
  fixed_matrix <- stats::model.matrix(design, data = model_frame)
  fixed_rank <- qr(fixed_matrix, tol = rank_tolerance)$rank
  coefficient_index <- if (is.character(coefficient)) {
    match(coefficient, colnames(fixed_matrix))
  } else {
    as.integer(coefficient)
  }
  if (length(coefficient_index) != 1L || is.na(coefficient_index) ||
      coefficient_index < 1L || coefficient_index > ncol(fixed_matrix)) {
    stop("coefficient is absent from the fixed-effect design.", call. = FALSE)
  }
  coefficient_name <- colnames(fixed_matrix)[coefficient_index]

  score <- prepared$score
  result <- data.frame(
    feature_id = colnames(score),
    coefficient = coefficient_name,
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    nonzero_units = colSums(score > 0),
    total_units = nrow(score),
    patient_count = length(unique(data$.patient_random_effect)),
    zero_fraction = colMeans(score <= 0),
    response_sd = NA_real_,
    singular = NA,
    response_scale = NA_real_,
    random_effect_variance = NA_real_,
    residual_variance = NA_real_,
    max_gradient = NA_real_,
    optimizer = NA_character_,
    optimizer_attempts = NA_character_,
    convergence_message = NA_character_,
    fit_error = NA_character_,
    status = "not_fitted",
    stringsAsFactors = FALSE
  )
  if (fixed_rank < ncol(fixed_matrix)) {
    result$status <- "rank_deficient"
  } else if (nrow(score) < min_units) {
    result$status <- "insufficient_units"
  } else {
    fixed_terms <- paste(deparse(design[[2L]]), collapse = "")
    model_formula <- stats::as.formula(paste(
      "response ~", fixed_terms, "+ (1 | .patient_random_effect)"
    ))
    for (feature in seq_len(ncol(score))) {
      if (result$nonzero_units[feature] < min_nonzero) {
        result$status[feature] <- "insufficient_nonzero"
        next
      }
      model_data <- data
      response <- log1p(score[, feature])
      response_mean <- mean(response)
      response_scale <- stats::sd(response)
      result$response_sd[feature] <- response_scale
      if (result$zero_fraction[feature] > max_zero_fraction) {
        result$status[feature] <- "insufficient_zero_support"
        next
      }
      if (!is.finite(response_scale) || response_scale <= min_response_sd) {
        result$status[feature] <- "insufficient_variation"
        next
      }
      model_data$response <- if (standardize_response) {
        (response - response_mean) / response_scale
      } else {
        response
      }
      fit_scale <- if (standardize_response) response_scale else 1
      result$response_scale[feature] <- fit_scale
      fit_with_optimizer <- function(name) {
        tryCatch(
          suppressWarnings(
            lme4::lmer(
              model_formula,
              data = model_data,
              REML = REML,
              control = lme4::lmerControl(
                optimizer = name,
                optCtrl = list(maxfun = 2e5),
                check.conv.singular = "ignore"
              )
            )
          ),
          error = function(error) error
        )
      }
      optimizer_sequence <- unique(c(optimizer, optimizer_fallbacks))
      errors <- character()
      attempted_optimizers <- character()
      selected_fit <- NULL
      selected_optimizer <- NA_character_
      selected_message <- NULL
      for (candidate_optimizer in optimizer_sequence) {
        attempted_optimizers <- c(attempted_optimizers, candidate_optimizer)
        candidate <- fit_with_optimizer(candidate_optimizer)
        if (inherits(candidate, "error")) {
          errors <- c(
            errors,
            paste0(candidate_optimizer, ": ", conditionMessage(candidate))
          )
          next
        }
        candidate_message <- candidate@optinfo$conv$lme4$messages
        if (!length(candidate_message)) {
          selected_fit <- candidate
          selected_optimizer <- candidate_optimizer
          selected_message <- NULL
          break
        }
        if (is.null(selected_fit)) {
          selected_fit <- candidate
          selected_optimizer <- candidate_optimizer
          selected_message <- candidate_message
        }
      }
      if (is.null(selected_fit)) {
        result$optimizer_attempts[feature] <- paste(
          attempted_optimizers, collapse = ","
        )
        result$fit_error[feature] <- paste(errors, collapse = "; ")
        numerical_error <- grepl(
          "Downdated VtV is not positive definite|not positive definite",
          result$fit_error[feature], ignore.case = TRUE
        )
        result$status[feature] <- if (numerical_error) {
          "numerical_failure"
        } else {
          "fit_failed"
        }
        next
      }
      fit <- selected_fit
      fit_optimizer <- selected_optimizer
      result$optimizer_attempts[feature] <- paste(
        attempted_optimizers, collapse = ","
      )
      coefficients <- tryCatch(
        summary(fit)$coefficients,
        error = function(error) NULL
      )
      if (is.null(coefficients) || !coefficient_name %in% rownames(coefficients) ||
          anyNA(coefficients[coefficient_name, seq_len(3L)])) {
        result$status[feature] <- "rank_deficient"
        next
      }
      estimate <- coefficients[coefficient_name, 1L] * fit_scale
      std_error <- coefficients[coefficient_name, 2L] * fit_scale
      statistic <- coefficients[coefficient_name, 3L]
      result$estimate[feature] <- estimate
      result$std_error[feature] <- std_error
      result$statistic[feature] <- statistic
      result$p_value[feature] <- 2 * stats::pnorm(-abs(statistic))
      result$singular[feature] <- lme4::isSingular(fit, tol = 1e-4)
      random_sd <- as.numeric(attr(lme4::VarCorr(fit)[[1L]], "stddev"))[1L]
      result$random_effect_variance[feature] <- (random_sd * fit_scale)^2
      result$residual_variance[feature] <- (stats::sigma(fit) * fit_scale)^2
      gradient <- fit@optinfo$derivs$gradient
      if (!is.null(gradient)) {
        result$max_gradient[feature] <- max(abs(gradient))
      }
      result$optimizer[feature] <- fit_optimizer
      convergence_message <- selected_message
      if (length(convergence_message)) {
        result$convergence_message[feature] <- paste(
          convergence_message, collapse = "; "
        )
        result$status[feature] <- "convergence_warning"
      } else {
        result$status[feature] <- "ok"
      }
    }
  }
  fdr_p_value <- result$p_value
  fdr_p_value[result$status != "ok"] <- NA_real_
  result$q_value <- stats::p.adjust(fdr_p_value, method = adjust_method)
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
