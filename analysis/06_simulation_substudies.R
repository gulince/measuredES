# 06_simulation_substudies.R
# Four focused sub-studies added in response to review:
#   A) finite-sample coverage of the rxx_se delta-method interval,
#      with alpha and its SE estimated per replicate (bootstrap SE)
#   B) unbalanced designs (B1) and heteroscedastic true scores (B2):
#      Hunter-Schmidt vs Bonett on their home turf
#   C) BCa bootstrap under markedly stronger skewness (chi-sq df = 1)
#   D) two-sided correction for correlations, both variables fallible

suppressMessages({
  library(measuredES)
  library(parallel)
})
source("analysis/00_helpers.R")

## ---------------------------------------------------------------- A
# rxx_se: does propagating the uncertainty of alpha-hat repair the
# undercoverage that short tests induce?

cells_a <- expand.grid(k_items = c(4L, 10L), n = c(30L, 60L),
                       rho = c(0.60, 0.85), stringsAsFactors = FALSE)
cells_a$cell_id <- seq_len(nrow(cells_a))
DELTA_A <- 0.50
R_A     <- 2000L
B_SE    <- 200L

se_alpha_boot <- function(X_mat, B) {
  n <- nrow(X_mat)
  a <- vapply(seq_len(B), function(b) {
    cronbach_alpha(X_mat[sample.int(n, n, replace = TRUE), , drop = FALSE])
  }, numeric(1))
  stats::sd(a, na.rm = TRUE)
}

run_a <- function(cell) {
  set.seed(20270000 + cell$cell_id)
  out <- matrix(NA_real_, R_A, 4,
                dimnames = list(NULL, c("cov_naive", "cov_se",
                                        "alpha_hat", "se_hat")))
  for (r in seq_len(R_A)) {
    dat <- generate_test_data(cell$n, cell$n, DELTA_A, cell$rho,
                              n_items = cell$k_items)
    X1 <- dat$X_total[dat$group == 1L]
    X2 <- dat$X_total[dat$group == 2L]
    a_hat <- max(min(cronbach_alpha(dat$X_mat), 0.999), 0.05)
    a_se  <- se_alpha_boot(dat$X_mat, B_SE)

    Sp    <- pooled_sd_fast(X1, X2)
    d_obs <- (mean(X2) - mean(X1)) / Sp
    ci0 <- measuredES:::.hs_ci_d(d_obs, a_hat, cell$n, cell$n, 0.95)
    ci1 <- measuredES:::.hs_ci_d(d_obs, a_hat, cell$n, cell$n, 0.95,
                                 rxx_se = a_se)
    out[r, ] <- c(ci0[1] <= DELTA_A && DELTA_A <= ci0[2],
                  ci1[1] <= DELTA_A && DELTA_A <= ci1[2],
                  a_hat, a_se)
  }
  data.frame(cell[, c("k_items", "n", "rho")],
             coverage_naive = mean(out[, "cov_naive"]),
             coverage_se    = mean(out[, "cov_se"]),
             mean_alpha     = mean(out[, "alpha_hat"]),
             mean_se_hat    = mean(out[, "se_hat"]),
             emp_sd_alpha   = stats::sd(out[, "alpha_hat"]))
}

## ---------------------------------------------------------------- B
# Unbalanced n (B1, equal variances) and heteroscedastic true scores
# (B2, balanced n). delta is defined on the average true-score SD,
# which equals the pooled/population SD in both blocks by design.

gen_het_target <- function(n1, n2, mean_diff, rho, sdT1 = 1, sdT2 = 1) {
  T1 <- rnorm(n1, 0,         sdT1)
  T2 <- rnorm(n2, mean_diff, sdT2)
  E1 <- rnorm(n1, 0, sdT1 * sqrt((1 - rho) / rho))
  E2 <- rnorm(n2, 0, sdT2 * sqrt((1 - rho) / rho))
  list(X1 = T1 + E1, X2 = T2 + E2)
}

cells_b <- rbind(
  expand.grid(block = "B1_unbalanced", n1 = c(30L, 60L), ratio_n = 3L,
              var_ratio = 1, pairing = "none", rho = c(0.60, 0.85),
              stringsAsFactors = FALSE),
  expand.grid(block = "B2_heteroscedastic", n1 = c(30L, 60L), ratio_n = 1L,
              var_ratio = 2, pairing = "none", rho = c(0.60, 0.85),
              stringsAsFactors = FALSE),
  expand.grid(block = "B3_combined", n1 = c(30L, 60L), ratio_n = 3L,
              var_ratio = 2, pairing = c("positive", "negative"),
              rho = c(0.60, 0.85), stringsAsFactors = FALSE)
)
cells_b$cell_id <- seq_len(nrow(cells_b))
DELTA_B <- 0.50
R_B     <- 5000L

# In B3, delta is defined on the population pooled SD (the parameter the
# pooled-d estimator targets), so coverage isolates interval behavior:
# positive pairing = larger variance in the LARGER group,
# negative pairing = larger variance in the SMALLER group.
run_b <- function(cell) {
  set.seed(20280000 + cell$cell_id)
  n1 <- cell$n1; n2 <- cell$n1 * cell$ratio_n
  if (cell$pairing == "negative") {
    sdT1 <- sqrt(cell$var_ratio); sdT2 <- 1   # big variance, small group
  } else {
    sdT1 <- 1; sdT2 <- sqrt(cell$var_ratio)
  }
  sd_pool <- sqrt(((n1 - 1) * sdT1^2 + (n2 - 1) * sdT2^2) /
                    (n1 + n2 - 2))
  out <- matrix(NA_real_, R_B, 4,
                dimnames = list(NULL, c("cov_hs", "wid_hs",
                                        "cov_bo", "wid_bo")))
  for (r in seq_len(R_B)) {
    dat <- gen_het_target(n1, n2, DELTA_B * sd_pool, cell$rho, sdT1, sdT2)
    X1 <- dat$X1; X2 <- dat$X2
    Sp    <- pooled_sd_fast(X1, X2)
    d_obs <- (mean(X2) - mean(X1)) / Sp
    ci_hs <- measuredES:::.hs_ci_d(d_obs, cell$rho, n2, n1, 0.95)
    ci_bo <- measuredES:::.bonett_ci_d(d_obs, cell$rho, n2, n1,
                                       stats::sd(X2), stats::sd(X1), 0.95)
    out[r, ] <- c(ci_hs[1] <= DELTA_B && DELTA_B <= ci_hs[2],
                  ci_hs[2] - ci_hs[1],
                  ci_bo[1] <= DELTA_B && DELTA_B <= ci_bo[2],
                  ci_bo[2] - ci_bo[1])
  }
  data.frame(cell[, c("block", "n1", "ratio_n", "var_ratio",
                      "pairing", "rho")],
             coverage_hs = mean(out[, "cov_hs"]),
             coverage_bo = mean(out[, "cov_bo"]),
             width_hs    = mean(out[, "wid_hs"]),
             width_bo    = mean(out[, "wid_bo"]))
}

## ---------------------------------------------------------------- C
# BCa under markedly stronger skewness: standardized chi-sq(1),
# skewness ~ 2.83 (vs ~ 1.41 in the main design).

gen_extreme <- function(n1, n2, delta, rho) {
  sf <- sqrt(2)   # sd of chi-sq(1)
  T1 <- (rchisq(n1, 1) - 1) / sf
  T2 <- (rchisq(n2, 1) - 1) / sf + delta
  sigma_E <- sqrt((1 - rho) / rho)
  E1 <- (rchisq(n1, 1) - 1) / sf * sigma_E
  E2 <- (rchisq(n2, 1) - 1) / sf * sigma_E
  list(X1 = T1 + E1, X2 = T2 + E2)
}

cells_c <- expand.grid(n = c(30L, 60L), rho = c(0.60, 0.85),
                       stringsAsFactors = FALSE)
cells_c$cell_id <- seq_len(nrow(cells_c))
DELTA_C <- 0.80
R_C     <- 1000L
B_C     <- 999L

run_c <- function(cell) {
  set.seed(20290000 + cell$cell_id)
  out <- matrix(NA_real_, R_C, 6,
                dimnames = list(NULL, c("cov_bca", "wid_bca", "cov_hs",
                                        "wid_hs", "cov_bo", "wid_bo")))
  for (r in seq_len(R_C)) {
    dat <- gen_extreme(cell$n, cell$n, DELTA_C, cell$rho)
    X1 <- dat$X1; X2 <- dat$X2
    ci_bca <- tryCatch(
      measuredES:::.boot_ci_d(X2, X1, cell$rho, 0.95, B_C),
      error = function(e) c(NA_real_, NA_real_))
    Sp    <- pooled_sd_fast(X1, X2)
    d_obs <- (mean(X2) - mean(X1)) / Sp
    ci_hs <- measuredES:::.hs_ci_d(d_obs, cell$rho, cell$n, cell$n, 0.95)
    ci_bo <- measuredES:::.bonett_ci_d(d_obs, cell$rho, cell$n, cell$n,
                                       stats::sd(X2), stats::sd(X1), 0.95)
    out[r, ] <- c(ci_bca[1] <= DELTA_C && DELTA_C <= ci_bca[2],
                  ci_bca[2] - ci_bca[1],
                  ci_hs[1] <= DELTA_C && DELTA_C <= ci_hs[2],
                  ci_hs[2] - ci_hs[1],
                  ci_bo[1] <= DELTA_C && DELTA_C <= ci_bo[2],
                  ci_bo[2] - ci_bo[1])
  }
  data.frame(cell[, c("n", "rho")],
             n_ok         = sum(!is.na(out[, "cov_bca"])),
             coverage_bca = mean(out[, "cov_bca"], na.rm = TRUE),
             coverage_hs  = mean(out[, "cov_hs"]),
             coverage_bo  = mean(out[, "cov_bo"]),
             width_bca    = mean(out[, "wid_bca"], na.rm = TRUE),
             width_hs     = mean(out[, "wid_hs"]),
             width_bo     = mean(out[, "wid_bo"]))
}

## ---------------------------------------------------------------- D
# Two-sided correction for correlations: both X and Y fallible.

cells_d <- expand.grid(rho_xy = c(0.30, 0.50), rel = c(0.70, 0.85),
                       n = c(100L, 250L), stringsAsFactors = FALSE)
cells_d$cell_id <- seq_len(nrow(cells_d))
R_D <- 5000L

run_d <- function(cell) {
  set.seed(20300000 + cell$cell_id)
  r_obs <- r_cor <- numeric(R_D)
  cov_r <- logical(R_D)
  sigE  <- sqrt((1 - cell$rel) / cell$rel)
  for (r in seq_len(R_D)) {
    Tx <- rnorm(cell$n)
    Ty <- cell$rho_xy * Tx + sqrt(1 - cell$rho_xy^2) * rnorm(cell$n)
    X  <- Tx + rnorm(cell$n, 0, sigE)
    Y  <- Ty + rnorm(cell$n, 0, sigE)
    r_obs[r] <- stats::cor(X, Y)
    fit <- mes_r(r_obs[r], n = cell$n, rxx = cell$rel, ryy = cell$rel)
    r_cor[r] <- fit$estimate
    cov_r[r] <- fit$ci[1] <= cell$rho_xy && cell$rho_xy <= fit$ci[2]
  }
  data.frame(cell[, c("rho_xy", "rel", "n")],
             bias_obs  = mean(r_obs) - cell$rho_xy,
             bias_corr = mean(r_cor) - cell$rho_xy,
             rmse_obs  = sqrt(mean((r_obs - cell$rho_xy)^2)),
             rmse_corr = sqrt(mean((r_cor - cell$rho_xy)^2)),
             coverage_fisher = mean(cov_r),
             pct_gt1   = 100 * mean(abs(r_cor) > 1))
}

## ------------------------------------------------------------ run
n_cores <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(cl, {
  suppressMessages(library(measuredES))
  source("analysis/00_helpers.R")
})
parallel::clusterExport(cl, c("DELTA_A", "R_A", "B_SE", "se_alpha_boot",
                              "DELTA_B", "R_B", "gen_het_target",
                              "DELTA_C", "R_C", "B_C", "gen_extreme",
                              "R_D"))

cat("Sub-studies started:", format(Sys.time(), "%H:%M:%S"), "\n")

res_a <- do.call(rbind, parallel::parLapply(
  cl, split(cells_a, cells_a$cell_id), run_a))
saveRDS(res_a, "analysis/results/substudy_rxxse.rds")
write.csv(res_a, "analysis/results/substudy_rxxse.csv", row.names = FALSE)
cat("A done:", format(Sys.time(), "%H:%M:%S"), "\n")
print(res_a, row.names = FALSE)

res_b <- do.call(rbind, parallel::parLapply(
  cl, split(cells_b, cells_b$cell_id), run_b))
saveRDS(res_b, "analysis/results/substudy_design.rds")
write.csv(res_b, "analysis/results/substudy_design.csv", row.names = FALSE)
cat("B done:", format(Sys.time(), "%H:%M:%S"), "\n")
print(res_b, row.names = FALSE)

res_c <- do.call(rbind, parallel::parLapply(
  cl, split(cells_c, cells_c$cell_id), run_c))
saveRDS(res_c, "analysis/results/substudy_extreme_skew.rds")
write.csv(res_c, "analysis/results/substudy_extreme_skew.csv",
          row.names = FALSE)
cat("C done:", format(Sys.time(), "%H:%M:%S"), "\n")
print(res_c, row.names = FALSE)

res_d <- do.call(rbind, parallel::parLapply(
  cl, split(cells_d, cells_d$cell_id), run_d))
saveRDS(res_d, "analysis/results/substudy_correlation.rds")
write.csv(res_d, "analysis/results/substudy_correlation.csv",
          row.names = FALSE)
cat("D done:", format(Sys.time(), "%H:%M:%S"), "\n")
print(res_d, row.names = FALSE)

parallel::stopCluster(cl)
cat("All sub-studies finished:", format(Sys.time(), "%H:%M:%S"), "\n")
