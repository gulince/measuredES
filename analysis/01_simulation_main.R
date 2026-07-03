# 01_simulation_main.R
# Study 1 (bias, RMSE) and Study 2 (CI coverage, width).
# Runs ~2-4 hours on 8 cores with R = 5000 replications per condition.

rm(list = ls())
set.seed(2026)

source("analysis/00_helpers.R")

suppressMessages({
  library(parallel)
  library(dplyr)
  library(tidyr)
})

n_cores <- max(1L, parallel::detectCores() - 1L)  # leave one core free
message(sprintf("Parallel cores available: %d", n_cores))


# ---- Factorial design ----
design <- expand.grid(
  delta             = c(0.20, 0.50, 0.80),
  n_per_group       = c(30L, 60L, 120L, 240L),
  rho_true          = c(0.60, 0.75, 0.85, 0.95),
  reliability_known = c(TRUE, FALSE),
  dist              = c("normal", "skewed"),
  stringsAsFactors  = FALSE
)
design$condition_id <- seq_len(nrow(design))

cat(sprintf("Total number of conditions: %d\n", nrow(design)))

R_REPS  <- 5000L      # replications per condition
N_ITEMS <- 10L        # number of items for alpha estimation


# ---- Run replications for one condition ----
run_condition <- function(cond_row, R = R_REPS, n_items = N_ITEMS) {

  reps <- vector("list", R)
  for (r in seq_len(R)) {
    out <- tryCatch(
      one_replication(n1 = cond_row$n_per_group,
                      n2 = cond_row$n_per_group,
                      delta = cond_row$delta,
                      rho_true = cond_row$rho_true,
                      reliability_known = cond_row$reliability_known,
                      dist = cond_row$dist,
                      n_items = n_items),
      error = function(e) NULL)
    if (!is.null(out)) reps[[r]] <- out
  }
  reps <- do.call(rbind, reps)

  delta_true <- cond_row$delta

  bias_obs  <- mean(reps$d_obs       - delta_true)
  bias_corr <- mean(reps$d_corrected - delta_true)

  rel_bias_obs  <- 100 * bias_obs  / delta_true
  rel_bias_corr <- 100 * bias_corr / delta_true

  rmse_obs  <- sqrt(mean((reps$d_obs       - delta_true)^2))
  rmse_corr <- sqrt(mean((reps$d_corrected - delta_true)^2))

  cov_hs <- mean(reps$ci_hs_low <= delta_true & delta_true <= reps$ci_hs_high)
  cov_bo <- mean(reps$ci_bo_low <= delta_true & delta_true <= reps$ci_bo_high)

  width_hs <- mean(reps$ci_hs_high - reps$ci_hs_low)
  width_bo <- mean(reps$ci_bo_high - reps$ci_bo_low)

  miss_low_hs  <- mean(reps$ci_hs_low  > delta_true)
  miss_high_hs <- mean(reps$ci_hs_high < delta_true)
  miss_low_bo  <- mean(reps$ci_bo_low  > delta_true)
  miss_high_bo <- mean(reps$ci_bo_high < delta_true)

  data.frame(
    condition_id      = cond_row$condition_id,
    delta             = cond_row$delta,
    n_per_group       = cond_row$n_per_group,
    rho_true          = cond_row$rho_true,
    reliability_known = cond_row$reliability_known,
    dist              = cond_row$dist,
    n_reps            = nrow(reps),

    # Study 1: point estimation
    bias_obs        = bias_obs,
    bias_corr       = bias_corr,
    rel_bias_obs    = rel_bias_obs,
    rel_bias_corr   = rel_bias_corr,
    rmse_obs        = rmse_obs,
    rmse_corr       = rmse_corr,
    mean_rho_hat    = mean(reps$rho_hat),

    # Study 2: CI performance
    coverage_hs     = cov_hs,
    coverage_bo     = cov_bo,
    width_hs        = width_hs,
    width_bo        = width_bo,
    miss_low_hs     = miss_low_hs,
    miss_high_hs    = miss_high_hs,
    miss_low_bo     = miss_low_bo,
    miss_high_bo    = miss_high_bo,

    stringsAsFactors = FALSE
  )
}


# ---- Run all conditions in parallel ----
cat("Starting simulation...\n")
tic <- Sys.time()

design_list <- split(design, design$condition_id)

# fork-based on unix, PSOCK cluster on Windows
if (n_cores > 1L && .Platform$OS.type == "unix") {
  results_list <- parallel::mclapply(design_list, run_condition,
                                     mc.cores = n_cores,
                                     mc.preschedule = FALSE)
} else {
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("R_REPS", "N_ITEMS", "one_replication",
                                "generate_two_groups", "generate_test_data",
                                "cronbach_alpha", "pooled_sd_fast",
                                "k_items_to_factor"))
  results_list <- parallel::parLapply(cl, design_list, run_condition)
  parallel::stopCluster(cl)
}

results <- do.call(rbind, results_list)
toc <- Sys.time()

cat(sprintf("Simulation complete. Elapsed: %s\n", format(toc - tic)))


# ---- Save ----
saveRDS(results, "analysis/results/simulation_main.rds")
write.csv(results, "analysis/results/simulation_main.csv",
          row.names = FALSE)

cat("Results saved to: analysis/results/simulation_main.{rds,csv}\n")
cat("Total rows: ", nrow(results), "\n")

cat("\n--- QUICK SUMMARY ---\n")
cat("Mean bias (uncorrected):",
    round(mean(results$bias_obs), 4), "\n")
cat("Mean bias (corrected):",
    round(mean(results$bias_corr), 4), "\n")
cat("Mean coverage (Hunter-Schmidt):",
    round(mean(results$coverage_hs), 3), "\n")
cat("Mean coverage (Bonett):",
    round(mean(results$coverage_bo), 3), "\n")
