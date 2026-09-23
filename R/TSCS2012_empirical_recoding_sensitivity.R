#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(foreign))

data_path <- Sys.getenv(
  "TSCS_DATA",
  file.path(getwd(), "data", "tscs2012q2.sav")
)
output_dir <- Sys.getenv(
  "TSCS_OUTPUT_DIR",
  file.path(getwd(), "outputs", "reviewer_simulation_results", "tscs2012_empirical")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) stop("TSCS data file not found: ", data_path)

raw <- read.spss(
  data_path,
  to.data.frame = TRUE,
  use.value.labels = FALSE,
  use.missings = FALSE,
  reencode = "CP950"
)

needed <- c("id", "a1", "a10", "b1", "j1", "j2", "j5", "k9", "wr")
missing_columns <- setdiff(needed, names(raw))
if (length(missing_columns)) {
  stop("Required variables missing: ", paste(missing_columns, collapse = ", "))
}

valid_j1 <- raw$j1 %in% 1:4
valid_j2 <- raw$j2 %in% 1:4
j1 <- ifelse(valid_j1, raw$j1, NA_real_)
j2 <- ifelse(valid_j2, raw$j2, NA_real_)
z_avg_raw <- rowMeans(cbind(j1, j2), na.rm = TRUE)
z_avg_raw[!valid_j1 & !valid_j2] <- NA_real_

# Original Read_TSCS2012(): reverse the average and compress half-step scores.
z_rev <- 5 - z_avg_raw
z4_original <- ifelse(
  is.na(z_rev), NA_real_,
  ifelse(z_rev >= 3.5, 4,
    ifelse(z_rev == 3.0, 3,
      ifelse(z_rev >= 2.0, 2, 1)
    )
  )
)
# Web Appendix C.2 uses different boundaries at raw averages 2.5 and 3.5.
z4_appendix <- ifelse(
  is.na(z_avg_raw), NA_real_,
  ifelse(z_avg_raw <= 1.5, 4,
    ifelse(z_avg_raw <= 2.5, 3,
      ifelse(z_avg_raw <= 3.5, 2, 1)
    )
  )
)
z3_original <- ifelse(
  is.na(z4_original), NA_real_,
  ifelse(z4_original == 4, 3, ifelse(z4_original == 3, 2, 1))
)
z3_appendix <- ifelse(
  is.na(z4_appendix), NA_real_,
  ifelse(z4_appendix == 4, 3, ifelse(z4_appendix == 3, 2, 1))
)
z2_original <- ifelse(
  is.na(z3_original), NA_real_, ifelse(z3_original == 1, 1, 2)
)

analysis_data <- data.frame(
  id = raw$id,
  rrt_yes = ifelse(raw$k9 == 1, 1, ifelse(raw$k9 == 2, 0, NA_real_)),
  male = ifelse(raw$a1 == 1, 1, ifelse(raw$a1 == 2, 0, NA_real_)),
  university = ifelse(raw$b1 %in% 16:21, 1,
    ifelse(raw$b1 %in% 1:15, 0, NA_real_)
  ),
  marital_group = ifelse(raw$a10 %in% 3:6, 1,
    ifelse(raw$a10 %in% 1:2, 0, NA_real_)
  ),
  j1 = j1,
  j2 = j2,
  j5 = ifelse(raw$j5 %in% 1:4, raw$j5, NA_real_),
  z_avg_raw = z_avg_raw,
  z4_original = z4_original,
  z3_original = z3_original,
  z2_original = z2_original,
  z4_appendix = z4_appendix,
  z3_appendix = z3_appendix,
  weight = suppressWarnings(as.numeric(raw$wr))
)

flow <- data.frame(
  step = c(
    "raw_records",
    "valid_K9_RRT_response",
    "valid_J1",
    "valid_J2",
    "valid_J5",
    "valid_gender",
    "valid_education",
    "valid_marital_group"
  ),
  n_remaining = NA_integer_,
  n_excluded_at_step = NA_integer_
)
keep <- rep(TRUE, nrow(analysis_data))
flow$n_remaining[1] <- sum(keep)
flow$n_excluded_at_step[1] <- 0L
rules <- list(
  !is.na(analysis_data$rrt_yes),
  !is.na(analysis_data$j1),
  !is.na(analysis_data$j2),
  !is.na(analysis_data$j5),
  !is.na(analysis_data$male),
  !is.na(analysis_data$university),
  !is.na(analysis_data$marital_group)
)
for (i in seq_along(rules)) {
  before <- sum(keep)
  keep <- keep & rules[[i]]
  flow$n_remaining[i + 1L] <- sum(keep)
  flow$n_excluded_at_step[i + 1L] <- before - sum(keep)
}
dat <- analysis_data[keep, ]

write.csv(flow, file.path(output_dir, "sample_exclusion_flow.csv"), row.names = FALSE)
write.csv(dat, file.path(output_dir, "analysis_data_reconstructed.csv"), row.names = FALSE)

freq_one <- function(x, coding, sample_name, sample_mask) {
  xx <- x[sample_mask]
  ww <- dat$weight[sample_mask]
  lev <- sort(unique(xx[!is.na(xx)]))
  total_w <- sum(ww[is.finite(ww)], na.rm = TRUE)
  do.call(rbind, lapply(lev, function(k) {
    in_cat <- xx == k
    data.frame(
      sample = sample_name,
      coding = coding,
      category = k,
      n = sum(in_cat, na.rm = TRUE),
      percent = 100 * mean(in_cat, na.rm = TRUE),
      weighted_n = sum(ww[in_cat & is.finite(ww)], na.rm = TRUE),
      weighted_percent = if (total_w > 0) {
        100 * sum(ww[in_cat & is.finite(ww)], na.rm = TRUE) / total_w
      } else NA_real_
    )
  }))
}

sample_masks <- list(
  full = rep(TRUE, nrow(dat)),
  ever_married = dat$marital_group == 1,
  non_married = dat$marital_group == 0
)
freqs <- list()
for (sample_name in names(sample_masks)) {
  for (coding in c(
    "j1_original", "j2_original", "z4_original_function",
    "z3_original_function", "z2_binary", "z4_web_appendix", "z3_web_appendix"
  )) {
    variable <- switch(coding,
      j1_original = dat$j1,
      j2_original = dat$j2,
      z4_original_function = dat$z4_original,
      z3_original_function = dat$z3_original,
      z2_binary = dat$z2_original,
      z4_web_appendix = dat$z4_appendix,
      z3_web_appendix = dat$z3_appendix
    )
    freqs[[length(freqs) + 1L]] <- freq_one(
      variable, coding, sample_name, sample_masks[[sample_name]]
    )
  }
}
frequency_table <- do.call(rbind, freqs)
write.csv(frequency_table, file.path(output_dir, "category_frequencies.csv"), row.names = FALSE)

pair_table <- as.data.frame(table(
  j1 = factor(dat$j1, levels = 1:4),
  j2 = factor(dat$j2, levels = 1:4),
  useNA = "ifany"
))
names(pair_table)[names(pair_table) == "Freq"] <- "n"
write.csv(pair_table, file.path(output_dir, "j1_by_j2_frequencies.csv"), row.names = FALSE)

inv_logit <- function(x) plogis(pmax(pmin(x, 35), -35))

decode_cutpoints <- function(block) {
  if (length(block) == 1L) return(block)
  c(block[1], block[1] + cumsum(exp(block[-1])))
}

encode_cutpoints <- function(cuts) {
  if (length(cuts) == 1L) return(cuts)
  c(cuts[1], log(pmax(diff(cuts), 1e-4)))
}

unpack_theta <- function(theta, K) {
  q <- K - 1L
  list(
    alpha = theta[1:3],
    cut0 = decode_cutpoints(theta[3 + seq_len(q)]),
    cut1 = decode_cutpoints(theta[3 + q + seq_len(q)]),
    gamma = theta[3 + 2 * q + 1:2]
  )
}

ordinal_matrix <- function(eta, cuts) {
  n <- length(eta)
  K <- length(cuts) + 1L
  # Match manuscript equation: logit Pr(Z <= k | Y, X) = beta_yk + X' gamma.
  cum <- sapply(cuts, function(cut) inv_logit(cut + eta))
  if (length(cuts) == 1L) cum <- matrix(cum, nrow = n, ncol = 1L)
  probs <- matrix(NA_real_, nrow = n, ncol = K)
  probs[, 1] <- cum[, 1]
  if (K > 2L) {
    for (k in 2:(K - 1L)) probs[, k] <- cum[, k] - cum[, k - 1L]
  }
  probs[, K] <- 1 - cum[, K - 1L]
  pmax(probs, 1e-12)
}

observed_components <- function(theta, z, r, X, K, p = 0.5, c = 0.25) {
  pars <- unpack_theta(theta, K)
  pi_y <- inv_logit(drop(X %*% pars$alpha))
  eta <- drop(X[, -1, drop = FALSE] %*% pars$gamma)
  pz0 <- ordinal_matrix(eta, pars$cut0)[cbind(seq_along(z), z)]
  pz1 <- ordinal_matrix(eta, pars$cut1)[cbind(seq_along(z), z)]
  kappa1 <- p + (1 - p) * c
  kappa0 <- (1 - p) * c
  pr0 <- ifelse(r == 1, kappa0, 1 - kappa0)
  pr1 <- ifelse(r == 1, kappa1, 1 - kappa1)
  comp0 <- (1 - pi_y) * pr0 * pz0
  comp1 <- pi_y * pr1 * pz1
  list(comp0 = comp0, comp1 = comp1, pi_y = pi_y)
}

nll <- function(theta, z, r, X, K) {
  comp <- observed_components(theta, z, r, X, K)
  value <- -sum(log(comp$comp0 + comp$comp1))
  if (!is.finite(value)) 1e100 else value
}

natural_vector <- function(theta, K) {
  pars <- unpack_theta(theta, K)
  out <- c(pars$alpha, pars$cut0, pars$cut1, pars$gamma)
  names(out) <- c(
    "alpha_intercept", "alpha_male", "alpha_university",
    paste0("beta0_cut", seq_len(K - 1L)),
    paste0("beta1_cut", seq_len(K - 1L)),
    "gamma_male", "gamma_university"
  )
  out
}

numeric_jacobian <- function(fun, x, eps = 1e-5) {
  y0 <- fun(x)
  J <- matrix(NA_real_, nrow = length(y0), ncol = length(x))
  for (j in seq_along(x)) {
    h <- eps * max(1, abs(x[j]))
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    J[, j] <- (fun(xp) - fun(xm)) / (2 * h)
  }
  rownames(J) <- names(y0)
  J
}

make_start <- function(z, r, X, K) {
  empirical_cuts <- qlogis(sapply(seq_len(K - 1L), function(k) mean(z <= k)))
  alpha0 <- qlogis(pmin(pmax((mean(r) - 0.125) / 0.5, 0.02), 0.98))
  c(
    alpha0, 0, 0,
    encode_cutpoints(empirical_cuts),
    encode_cutpoints(empirical_cuts + 0.15),
    0, 0
  )
}

fit_joint <- function(d, z_name, sample_name, n_starts = 12L) {
  z <- as.integer(d[[z_name]])
  K <- max(z)
  r <- as.integer(d$rrt_yes)
  X <- cbind(1, d$male, d$university)
  start <- make_start(z, r, X, K)
  set.seed(20260920 + K + nrow(d))
  fits <- vector("list", n_starts)
  for (s in seq_len(n_starts)) {
    initial <- if (s == 1L) start else start + rnorm(length(start), 0, 0.35)
    fits[[s]] <- tryCatch(
      optim(
        initial, nll, z = z, r = r, X = X, K = K,
        method = "BFGS", hessian = TRUE,
        control = list(maxit = 3000, reltol = 1e-11)
      ),
      error = function(e) NULL
    )
  }
  valid <- which(vapply(fits, function(x) {
    !is.null(x) && is.finite(x$value) && all(is.finite(x$par))
  }, logical(1)))
  if (!length(valid)) stop("All starts failed for ", sample_name, " / ", z_name)
  best_index <- valid[which.min(vapply(fits[valid], `[[`, numeric(1), "value"))]
  fit <- fits[[best_index]]
  natural <- natural_vector(fit$par, K)

  eig <- tryCatch(eigen(fit$hessian, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, length(fit$par))
  )
  positive_hessian <- all(is.finite(eig)) && min(eig) > 1e-7
  se <- rep(NA_real_, length(natural))
  if (positive_hessian) {
    covariance <- tryCatch(solve(fit$hessian), error = function(e) NULL)
    if (!is.null(covariance)) {
      J <- numeric_jacobian(function(x) natural_vector(x, K), fit$par)
      vnat <- J %*% covariance %*% t(J)
      se <- sqrt(pmax(diag(vnat), 0))
    }
  }
  names(se) <- names(natural)

  close_values <- valid[vapply(fits[valid], function(x) abs(x$value - fit$value) < 1e-5, logical(1))]
  max_solution_difference <- if (length(close_values) > 1L) {
    max(vapply(fits[close_values], function(x) max(abs(x$par - fit$par)), numeric(1)))
  } else 0

  comp <- observed_components(fit$par, z, r, X, K)
  posterior <- comp$comp1 / (comp$comp0 + comp$comp1)
  summary_row <- data.frame(
    sample = sample_name,
    coding = z_name,
    categories = K,
    n = nrow(d),
    rrt_yes_n = sum(r == 1),
    converged = fit$convergence == 0,
    n_starts = n_starts,
    finite_starts = length(valid),
    starts_at_best_likelihood = length(close_values),
    max_parameter_difference_at_best = max_solution_difference,
    positive_definite_hessian = positive_hessian,
    log_likelihood = -fit$value,
    n_parameters = length(fit$par),
    AIC = 2 * length(fit$par) + 2 * fit$value,
    BIC = log(nrow(d)) * length(fit$par) + 2 * fit$value,
    marginal_prevalence = mean(comp$pi_y),
    posterior_mean_prevalence = mean(posterior),
    posterior_sd = sd(posterior)
  )
  parameter_rows <- data.frame(
    sample = sample_name,
    coding = z_name,
    parameter = names(natural),
    estimate = as.numeric(natural),
    SE = as.numeric(se),
    z_value = as.numeric(natural / se),
    lower_95 = as.numeric(natural - 1.96 * se),
    upper_95 = as.numeric(natural + 1.96 * se)
  )
  posterior_rows <- do.call(rbind, lapply(sort(unique(z)), function(k) {
    take <- z == k
    data.frame(
      sample = sample_name,
      coding = z_name,
      category = k,
      n = sum(take),
      mean_posterior_EMS = mean(posterior[take]),
      sd_posterior_EMS = sd(posterior[take]),
      mean_model_prevalence = mean(comp$pi_y[take])
    )
  }))
  list(summary = summary_row, parameters = parameter_rows, posterior = posterior_rows)
}

fit_results <- list()
for (sample_name in names(sample_masks)) {
  dsub <- dat[sample_masks[[sample_name]], ]
  for (coding in c(
    "z3_original", "z4_original", "z2_original",
    "z3_appendix", "z4_appendix"
  )) {
    message("Fitting ", sample_name, " / ", coding, " (n=", nrow(dsub), ")")
    fit_results[[paste(sample_name, coding, sep = "_")]] <- fit_joint(
      dsub, coding, sample_name
    )
  }
}

model_summary <- do.call(rbind, lapply(fit_results, `[[`, "summary"))
parameter_estimates <- do.call(rbind, lapply(fit_results, `[[`, "parameters"))
posterior_by_category <- do.call(rbind, lapply(fit_results, `[[`, "posterior"))

write.csv(model_summary, file.path(output_dir, "model_summary.csv"), row.names = FALSE)
write.csv(parameter_estimates, file.path(output_dir, "parameter_estimates.csv"), row.names = FALSE)
write.csv(posterior_by_category, file.path(output_dir, "posterior_EMS_by_attitude_category.csv"), row.names = FALSE)

audit <- data.frame(
  item = c(
    "source_file", "raw_n", "reconstructed_analysis_n",
    "reconstructed_ever_married_n", "reconstructed_non_married_n",
    "manuscript_full_n", "manuscript_ever_married_n", "manuscript_non_married_n",
    "difference_full_n", "p", "c", "education_rule", "marital_rule"
  ),
  value = c(
    data_path, nrow(raw), nrow(dat),
    sum(dat$marital_group == 1), sum(dat$marital_group == 0),
    1838, 1236, 602, nrow(dat) - 1838,
    0.5, 0.25, "b1 in 16:21", "a10 in 3:6 vs a10 in 1:2"
  )
)
write.csv(audit, file.path(output_dir, "data_audit.csv"), row.names = FALSE)

saveRDS(
  list(
    audit = audit,
    exclusion_flow = flow,
    frequencies = frequency_table,
    model_summary = model_summary,
    parameters = parameter_estimates,
    posterior = posterior_by_category
  ),
  file.path(output_dir, "tscs2012_empirical_sensitivity_results.rds")
)

message("Completed. Results written to: ", output_dir)
