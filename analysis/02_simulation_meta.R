# 02_simulation_meta.R
# Study 3: between-study reliability heterogeneity in meta-analytic
# point estimation. Each of the k synthesized studies has its own
# reliability and sample size.

rm(list = ls())
# Reproducibility: seeds are set per condition inside run_meta_condition()
# (base 20262000 + condition_id); results do not depend on core count.

source("analysis/00_helpers.R")

suppressMessages({
  library(parallel)
  library(metafor)
})

n_cores <- as.integer(Sys.getenv("MES_CORES",
                                 max(1L, parallel::detectCores() - 1L)))


# ---- Design ----
design_meta <- expand.grid(
  k_studies     = c(5L, 10L, 20L),
  delta_mean    = c(0.30, 0.50),
  rho_pattern   = c("homogeneous", "heterogeneous"),
  stringsAsFactors = FALSE
)
design_meta$condition_id <- seq_len(nrow(design_meta))

R_REPS_META <- 1000L     # replications per condition


simulate_one_meta <- function(k, delta_mean, rho_pattern) {

  n_per_study <- sample(c(40L, 60L, 100L, 150L, 250L), k, replace = TRUE)

  rho_per_study <- if (rho_pattern == "homogeneous") {
    rep(0.85, k)
  } else {
    # Beta(8, 4): mean ~ 0.67, sd ~ 0.13, clipped to [0.40, 0.95]
    pmin(pmax(stats::rbeta(k, 8, 4), 0.40), 0.95)
  }

  # random-effects: small heterogeneity in the true deltas
  delta_per_study <- stats::rnorm(k, mean = delta_mean, sd = 0.10)

  yi_obs    <- numeric(k)
  vi_obs    <- numeric(k)
  yi_corr   <- numeric(k)
  vi_corr   <- numeric(k)
  rho_known <- rho_per_study   # Study 3 assumes reliability is known

  for (i in seq_len(k)) {
    dat <- generate_two_groups(n_per_study[i], n_per_study[i],
                               delta_per_study[i], rho_per_study[i],
                               dist = "normal")
    Sp <- pooled_sd_fast(dat$X1, dat$X2)
    if (Sp == 0) return(NULL)

    d_obs <- (mean(dat$X2) - mean(dat$X1)) / Sp
    n1 <- length(dat$X1); n2 <- length(dat$X2)

    var_d_obs <- (n1 + n2) / (n1 * n2) + d_obs^2 / (2 * (n1 + n2))

    yi_obs[i]  <- d_obs
    vi_obs[i]  <- var_d_obs
    yi_corr[i] <- d_obs / sqrt(rho_known[i])
    vi_corr[i] <- var_d_obs / rho_known[i]
  }

  fit_obs  <- tryCatch(metafor::rma(yi = yi_obs, vi = vi_obs,
                                    method = "REML"),
                       error = function(e) NULL)
  fit_corr <- tryCatch(metafor::rma(yi = yi_corr, vi = vi_corr,
                                    method = "REML"),
                       error = function(e) NULL)

  if (is.null(fit_obs) || is.null(fit_corr)) return(NULL)

  data.frame(
    delta_true        = delta_mean,
    pooled_obs        = as.numeric(fit_obs$b),
    pooled_corr       = as.numeric(fit_corr$b),
    se_obs            = fit_obs$se,
    se_corr           = fit_corr$se,
    tau2_obs          = fit_obs$tau2,
    tau2_corr         = fit_corr$tau2,
    ci_obs_low        = fit_obs$ci.lb,
    ci_obs_high       = fit_obs$ci.ub,
    ci_corr_low       = fit_corr$ci.lb,
    ci_corr_high      = fit_corr$ci.ub
  )
}


run_meta_condition <- function(cond_row, R = R_REPS_META) {

  set.seed(20262000L + cond_row$condition_id)
  reps <- vector("list", R)
  for (r in seq_len(R)) {
    reps[[r]] <- simulate_one_meta(cond_row$k_studies,
                                   cond_row$delta_mean,
                                   cond_row$rho_pattern)
  }
  reps <- do.call(rbind, reps[!sapply(reps, is.null)])

  delta_true <- cond_row$delta_mean

  data.frame(
    condition_id  = cond_row$condition_id,
    k_studies     = cond_row$k_studies,
    delta_mean    = delta_true,
    rho_pattern   = cond_row$rho_pattern,
    n_reps        = nrow(reps),

    bias_pooled_obs   = mean(reps$pooled_obs)   - delta_true,
    bias_pooled_corr  = mean(reps$pooled_corr)  - delta_true,
    rmse_pooled_obs   = sqrt(mean((reps$pooled_obs  - delta_true)^2)),
    rmse_pooled_corr  = sqrt(mean((reps$pooled_corr - delta_true)^2)),

    coverage_obs   = mean(reps$ci_obs_low  <= delta_true &
                          delta_true <= reps$ci_obs_high),
    coverage_corr  = mean(reps$ci_corr_low <= delta_true &
                          delta_true <= reps$ci_corr_high),

    mean_tau2_obs   = mean(reps$tau2_obs),
    mean_tau2_corr  = mean(reps$tau2_corr),

    stringsAsFactors = FALSE
  )
}


# ---- Run ----
cat("Starting Study 3 (meta-analytic simulation)...\n")
tic <- Sys.time()

design_meta_list <- split(design_meta, design_meta$condition_id)

if (n_cores > 1L && .Platform$OS.type == "unix") {
  results_meta <- parallel::mclapply(design_meta_list, run_meta_condition,
                                     mc.cores = n_cores,
                                     mc.preschedule = FALSE)
} else {
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("R_REPS_META", "simulate_one_meta",
                                "generate_two_groups", "pooled_sd_fast"))
  parallel::clusterEvalQ(cl, library(metafor))
  results_meta <- parallel::parLapply(cl, design_meta_list, run_meta_condition)
  parallel::stopCluster(cl)
}

results_meta <- do.call(rbind, results_meta)
toc <- Sys.time()

cat(sprintf("Complete. Elapsed: %s\n", format(toc - tic)))

saveRDS(results_meta, "analysis/results/simulation_meta.rds")
write.csv(results_meta, "analysis/results/simulation_meta.csv",
          row.names = FALSE)

cat("Results saved to: analysis/results/simulation_meta.{rds,csv}\n")
print(results_meta)
