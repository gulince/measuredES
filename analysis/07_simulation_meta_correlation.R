# 07_simulation_meta_correlation.R
# Sub-study E: meta-analytic pooling of corrected correlations. Mirrors
# Study 3 (02_simulation_meta.R) but for mes_r()/mes_meta() instead of
# mes_d(), directly exercising the mes_meta() correlation variance
#   Var(r*) = (1 - r^2)^2 / ((n - 1) * rxx * ryy)
# added in response to review, in an actual metafor::rma() pooling step.

rm(list = ls())
set.seed(2026)

if (!requireNamespace("measuredES", quietly = TRUE)) {
  devtools::load_all(".")
} else {
  library(measuredES)
}
suppressMessages({
  library(parallel)
  library(metafor)
})

n_cores <- max(1L, parallel::detectCores() - 1L)


# ---- Design (mirrors Study 3's k/reliability-pattern crossing) ----
design_meta_r <- expand.grid(
  k_studies   = c(5L, 10L, 20L),
  rho_xy_mean = c(0.30, 0.50),
  rel_pattern = c("homogeneous", "heterogeneous"),
  stringsAsFactors = FALSE
)
design_meta_r$condition_id <- seq_len(nrow(design_meta_r))

R_REPS_META_R <- 1000L    # replications per condition


simulate_one_meta_r <- function(k, rho_xy_mean, rel_pattern) {

  n_per_study <- sample(c(40L, 60L, 100L, 150L, 250L), k, replace = TRUE)

  rel_per_study <- if (rel_pattern == "homogeneous") {
    rep(0.85, k)
  } else {
    # Beta(8, 4): mean ~ 0.67, sd ~ 0.13, clipped -- identical to Study 3
    pmin(pmax(stats::rbeta(k, 8, 4), 0.40), 0.95)
  }

  # random-effects: small between-study heterogeneity in the true r,
  # clipped well short of +-1 so the bivariate generator stays well-defined
  rho_xy_per_study <- pmin(pmax(
    stats::rnorm(k, mean = rho_xy_mean, sd = 0.10), -0.90), 0.90)

  r_obs  <- numeric(k)
  fits_r <- vector("list", k)

  for (i in seq_len(k)) {
    n_i <- n_per_study[i]
    Tx  <- stats::rnorm(n_i)
    Ty  <- rho_xy_per_study[i] * Tx +
      sqrt(1 - rho_xy_per_study[i]^2) * stats::rnorm(n_i)

    sigma_e <- sqrt((1 - rel_per_study[i]) / rel_per_study[i])
    X <- Tx + stats::rnorm(n_i, 0, sigma_e)
    Y <- Ty + stats::rnorm(n_i, 0, sigma_e)

    r_obs[i]    <- stats::cor(X, Y)
    fits_r[[i]] <- mes_r(r_obs[i], n = n_i,
                         rxx = rel_per_study[i], ryy = rel_per_study[i],
                         ci = FALSE)
  }

  # the actual shipped mes_meta() code path -- this is what the sub-study
  # is validating, not a re-implementation of its formula
  meta_df <- mes_meta(fits_r)

  # uncorrected baseline: same large-sample variance, no reliability scaling
  vi_obs <- (1 - r_obs^2)^2 / (n_per_study - 1)

  fit_obs  <- tryCatch(
    metafor::rma(yi = r_obs, vi = vi_obs, method = "REML"),
    error = function(e) NULL)
  fit_corr <- tryCatch(
    metafor::rma(yi = meta_df$yi, vi = meta_df$vi, method = "REML"),
    error = function(e) NULL)

  if (is.null(fit_obs) || is.null(fit_corr)) return(NULL)

  data.frame(
    rho_xy_true  = rho_xy_mean,
    pooled_obs   = as.numeric(fit_obs$b),
    pooled_corr  = as.numeric(fit_corr$b),
    ci_obs_low   = fit_obs$ci.lb,
    ci_obs_high  = fit_obs$ci.ub,
    ci_corr_low  = fit_corr$ci.lb,
    ci_corr_high = fit_corr$ci.ub
  )
}


run_meta_r_condition <- function(cond_row, R = R_REPS_META_R) {
  reps <- vector("list", R)
  for (r in seq_len(R)) {
    reps[[r]] <- simulate_one_meta_r(cond_row$k_studies,
                                     cond_row$rho_xy_mean,
                                     cond_row$rel_pattern)
  }
  reps <- do.call(rbind, reps[!vapply(reps, is.null, logical(1L))])

  rho_true <- cond_row$rho_xy_mean

  data.frame(
    condition_id = cond_row$condition_id,
    k_studies    = cond_row$k_studies,
    rho_xy_mean  = rho_true,
    rel_pattern  = cond_row$rel_pattern,
    n_reps       = nrow(reps),

    bias_pooled_obs  = mean(reps$pooled_obs)  - rho_true,
    bias_pooled_corr = mean(reps$pooled_corr) - rho_true,
    rmse_pooled_obs  = sqrt(mean((reps$pooled_obs  - rho_true)^2)),
    rmse_pooled_corr = sqrt(mean((reps$pooled_corr - rho_true)^2)),

    coverage_obs  = mean(reps$ci_obs_low  <= rho_true &
                         rho_true <= reps$ci_obs_high),
    coverage_corr = mean(reps$ci_corr_low <= rho_true &
                         rho_true <= reps$ci_corr_high)
  )
}


# ---- Run ----
cat("Starting correlation meta-analysis sub-study...\n")
tic <- Sys.time()

design_list <- split(design_meta_r, design_meta_r$condition_id)

cl <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(cl, {
  suppressMessages(library(measuredES))
  suppressMessages(library(metafor))
})
parallel::clusterExport(cl, c("R_REPS_META_R", "simulate_one_meta_r"))
results_meta_r <- parallel::parLapply(cl, design_list, run_meta_r_condition)
parallel::stopCluster(cl)

results_meta_r <- do.call(rbind, results_meta_r)
toc <- Sys.time()
cat(sprintf("Complete. Elapsed: %s\n", format(toc - tic)))

saveRDS(results_meta_r, "analysis/results/substudy_meta_correlation.rds")
write.csv(results_meta_r, "analysis/results/substudy_meta_correlation.csv",
          row.names = FALSE)

cat("Results saved to: analysis/results/substudy_meta_correlation.{rds,csv}\n")
print(results_meta_r, row.names = FALSE)
