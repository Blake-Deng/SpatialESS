# Experimental cohort robustness layer for SpatialESS 0.1.4.
# This module consumes saved communication scores; it does not change inference.
audit_cohort_designs <- function(metadata, designs) {
  do.call(rbind, lapply(names(designs), function(name) {
    matrix <- model.matrix(designs[[name]], metadata)
    rank <- qr(matrix, tol = 1e-10)$rank
    data.frame(design = name, rank = rank, columns = ncol(matrix),
               estimable = rank == ncol(matrix))
  }))
}

fit_patient_run_hc3 <- function(prepared, control = "Control", disease = "Disease") {
  metadata <- prepared$unit_metadata
  required <- c("patient_id", "condition", "run")
  stopifnot(all(required %in% names(metadata)), !anyNA(metadata[, required]),
            !anyNA(prepared$score), all(prepared$score >= 0),
            nrow(prepared$score) == nrow(metadata))
  patients <- unique(as.character(metadata$patient_id))
  runs <- sort(unique(as.character(metadata$run)))
  response <- log1p(prepared$score)
  y <- t(vapply(patients, function(patient) {
    selected <- metadata$patient_id == patient
    stopifnot(length(unique(metadata$condition[selected])) == 1L)
    colMeans(response[selected, , drop = FALSE])
  }, numeric(ncol(response))))
  rownames(y) <- patients; colnames(y) <- colnames(response)
  condition <- vapply(patients, function(p) as.character(metadata$condition[metadata$patient_id == p][1]), character(1))
  stopifnot(setequal(unique(condition), c(control, disease)))
  x <- cbind(intercept = 1, condition = as.numeric(condition == disease))
  if (length(runs) > 1L) {
    fractions <- vapply(runs[-1], function(run) {
      vapply(patients, function(p) mean(metadata$run[metadata$patient_id == p] == run), numeric(1))
    }, numeric(length(patients)))
    x <- cbind(x, fractions)
  }
  if (qr(x, tol = 1e-10)$rank < ncol(x)) stop("Patient run-composition design is not identifiable; no automatic covariate removal.")
  df <- nrow(x) - ncol(x)
  if (df < 3L) stop("Insufficient residual degrees of freedom for patient HC3 sensitivity.")
  inverse <- solve(crossprod(x)); weights <- inverse %*% t(x)
  coefficients <- weights %*% y
  leverage <- rowSums(x * (x %*% inverse))
  if (any(leverage >= 1 - 1e-10)) stop("Unit leverage prevents HC3 estimation.")
  residual <- y - x %*% coefficients
  variance <- as.numeric((weights[2, ]^2 / (1-leverage)^2) %*% (residual^2))
  se <- sqrt(pmax(variance, 0))
  nonzero <- colSums(y > 0)
  eligible <- nonzero >= 2L & apply(y, 2, sd) > 1e-12 & is.finite(se) & se > 0
  p <- rep(NA_real_, ncol(y))
  p[eligible] <- 2 * pt(-abs(coefficients[2, eligible] / se[eligible]), df)
  result <- data.frame(feature_id = colnames(y), estimate = coefficients[2, ],
    std_error = se, ci_lower = coefficients[2, ] - qt(.975, df)*se,
    ci_upper = coefficients[2, ] + qt(.975, df)*se, p_value = p,
    q_value = p.adjust(p, "BH"), residual_df = df, total_patients = nrow(y),
    nonzero_patients = nonzero,
    status = ifelse(eligible, "ok", "insufficient_support_or_variation"), row.names = NULL)
  list(results = result, patient_response = y,
       patient_design = data.frame(patient_id = patients, condition = condition, x, row.names = NULL),
       parameters = list(response = "mean(log1p(slice_score)), including zeros",
                         se = "HC3", test = "t with n_patient - rank(X) degrees of freedom",
                         run_reference = runs[1], adjustment = "run composition fractions",
                         fdr = "BH across all finite eligible p values", role = "sensitivity"))
}

classify_cohort_evidence <- function(primary, run_adjusted, patient_hc3, alpha = .05) {
  # These labels describe agreement, not a new calibrated test or biological truth.
  keys <- primary$feature_id
  align <- function(table) {
    stopifnot(!anyDuplicated(table$feature_id), setequal(keys, table$feature_id))
    table[match(keys, table$feature_id), , drop = FALSE]
  }
  stopifnot(!anyDuplicated(keys))
  run_adjusted <- align(run_adjusted); patient_hc3 <- align(patient_hc3)
  available <- function(x) x$status == "ok" & is.finite(x$q_value) & is.finite(x$estimate)
  is_hit <- function(x) available(x) & x$q_value < alpha
  p_hit <- is_hit(primary); run_hit <- is_hit(run_adjusted); patient_hit <- is_hit(patient_hc3)
  direction <- sign(primary$estimate) == sign(run_adjusted$estimate) &
               sign(primary$estimate) == sign(patient_hc3$estimate)
  direction[is.na(direction)] <- FALSE
  all_available <- available(primary) & available(run_adjusted) & available(patient_hc3)
  label <- rep("not_primary_hit", length(keys))
  label[!available(primary)] <- "primary_fit_unavailable"
  label[p_hit & !all_available] <- "primary_hit_sensitivity_unavailable"
  label[p_hit & all_available] <- "primary_hit_sensitivity_dependent"
  label[p_hit & run_hit & patient_hit & direction] <- "supported_by_run_and_patient_sensitivity"
  data.frame(feature_id = keys, primary_q = primary$q_value, run_adjusted_q = run_adjusted$q_value,
    patient_hc3_q = patient_hc3$q_value, all_models_available = all_available,
    same_direction = direction, evidence_label = label,
    primary_singular = primary$singular, row.names = NULL)
}
