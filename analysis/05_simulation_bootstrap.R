# 05_simulation_bootstrap.R
# Focused sub-study: coverage and width of the BCa bootstrap interval,
# compared with Hunter-Schmidt and Bonett on the same replicates.
# Cells target the setting where BCa should matter most: skewed scores
# and small samples (with normal cells as the reference).

suppressMessages({
  library(measuredES)
  library(parallel)
})
source("analysis/00_helpers.R")

cells <- expand.grid(
  dist  = c("normal", "skewed"),
  n     = c(30L, 60L),
  rho   = c(0.60, 0.85),
  delta = c(0.50, 0.80),
  stringsAsFactors = FALSE
)
cells$cell_id <- seq_len(nrow(cells))

R <- 1000L
B <- 999L

run_cell <- function(cell) {
  set.seed(20260000 + cell$cell_id)
  DELTA <- cell$delta
  out <- matrix(NA_real_, nrow = R, ncol = 6,
                dimnames = list(NULL, c("cov_bca", "wid_bca",
                                        "cov_hs",  "wid_hs",
                                        "cov_bo",  "wid_bo")))
  for (r in seq_len(R)) {
    dat <- generate_two_groups(cell$n, cell$n, DELTA, cell$rho, cell$dist)
    X1 <- dat$X1; X2 <- dat$X2   # delta lives on X2

    ci_bca <- tryCatch(
      measuredES:::.boot_ci_d(X2, X1, cell$rho, 0.95, B),
      error = function(e) c(NA_real_, NA_real_))

    Sp    <- pooled_sd_fast(X1, X2)
    d_obs <- (mean(X2) - mean(X1)) / Sp
    ci_hs <- measuredES:::.hs_ci_d(d_obs, cell$rho, cell$n, cell$n, 0.95)
    ci_bo <- measuredES:::.bonett_ci_d(d_obs, cell$rho, cell$n, cell$n,
                                       stats::sd(X2), stats::sd(X1), 0.95)

    out[r, ] <- c(ci_bca[1] <= DELTA && DELTA <= ci_bca[2],
                  ci_bca[2] - ci_bca[1],
                  ci_hs[1]  <= DELTA && DELTA <= ci_hs[2],
                  ci_hs[2]  - ci_hs[1],
                  ci_bo[1]  <= DELTA && DELTA <= ci_bo[2],
                  ci_bo[2]  - ci_bo[1])
  }
  data.frame(cell[, c("dist", "n", "rho", "delta")],
             n_ok         = sum(!is.na(out[, "cov_bca"])),
             coverage_bca = mean(out[, "cov_bca"], na.rm = TRUE),
             coverage_hs  = mean(out[, "cov_hs"],  na.rm = TRUE),
             coverage_bo  = mean(out[, "cov_bo"],  na.rm = TRUE),
             width_bca    = mean(out[, "wid_bca"], na.rm = TRUE),
             width_hs     = mean(out[, "wid_hs"],  na.rm = TRUE),
             width_bo     = mean(out[, "wid_bo"],  na.rm = TRUE))
}

n_cores <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(min(n_cores, nrow(cells)))
parallel::clusterEvalQ(cl, {
  suppressMessages(library(measuredES))
  source("analysis/00_helpers.R")
})
parallel::clusterExport(cl, c("R", "B", "run_cell"))

cat("Bootstrap sub-study started:", format(Sys.time(), "%H:%M:%S"), "\n")
res <- parallel::parLapply(cl, split(cells, cells$cell_id), run_cell)
parallel::stopCluster(cl)

results_boot <- do.call(rbind, res)
saveRDS(results_boot, "analysis/results/simulation_bootstrap.rds")
write.csv(results_boot, "analysis/results/simulation_bootstrap.csv",
          row.names = FALSE)
cat("Done:", format(Sys.time(), "%H:%M:%S"), "\n")
print(results_boot, row.names = FALSE)
