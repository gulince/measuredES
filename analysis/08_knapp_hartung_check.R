# 08_knapp_hartung_check.R
# Supplementary check for Section 5.5 of the manuscript: does the Knapp-Hartung
# adjustment (test = "knha" in metafor::rma) close the coverage gap of the
# pooled corrected correlation? Uses the same data-generating process as
# 07_simulation_meta_correlation.R and compares the Wald and Knapp-Hartung
# intervals computed on the SAME simulated replicates, for all 12 conditions.
#
# Run from the package root:  source("analysis/08_knapp_hartung_check.R")
# Output: analysis/results/substudy_knapp_hartung.{rds,csv}

rm(list = ls())
set.seed(20261011)

if (!requireNamespace("measuredES", quietly = TRUE)) {
  devtools::load_all(".")
} else {
  library(measuredES)
}
suppressMessages(library(metafor))

R_REPS <- 1000L

design <- expand.grid(
  k_studies   = c(5L, 10L, 20L),
  rho_xy_mean = c(0.30, 0.50),
  rel_pattern = c("homogeneous", "heterogeneous"),
  stringsAsFactors = FALSE
)

# One replicate: same generator as 07_simulation_meta_correlation.R
one_rep <- function(k, rho_xy_mean, rel_pattern) {
  n_per_study <- sample(c(40L, 60L, 100L, 150L, 250L), k, replace = TRUE)
  rel <- if (rel_pattern == "homogeneous") rep(0.85, k) else
    pmin(pmax(stats::rbeta(k, 8, 4), 0.40), 0.95)
  rho_i <- pmin(pmax(stats::rnorm(k, rho_xy_mean, 0.10), -0.90), 0.90)

  fits <- vector("list", k)
  for (i in seq_len(k)) {
    n_i <- n_per_study[i]
    Tx <- stats::rnorm(n_i)
    Ty <- rho_i[i] * Tx + sqrt(1 - rho_i[i]^2) * stats::rnorm(n_i)
    se <- sqrt((1 - rel[i]) / rel[i])
    X <- Tx + stats::rnorm(n_i, 0, se)
    Y <- Ty + stats::rnorm(n_i, 0, se)
    fits[[i]] <- mes_r(stats::cor(X, Y), n = n_i, rxx = rel[i], ryy = rel[i],
                       ci = FALSE)
  }
  md <- mes_meta(fits)

  wald <- tryCatch(metafor::rma(yi = md$yi, vi = md$vi, method = "REML"),
                   error = function(e) NULL)
  knha <- tryCatch(metafor::rma(yi = md$yi, vi = md$vi, method = "REML",
                                test = "knha"),
                   error = function(e) NULL)
  if (is.null(wald) || is.null(knha)) return(NULL)
  c(wald_lo = wald$ci.lb, wald_hi = wald$ci.ub,
    knha_lo = knha$ci.lb, knha_hi = knha$ci.ub)
}

out <- lapply(seq_len(nrow(design)), function(j) {
  d <- design[j, ]
  reps <- do.call(rbind, lapply(seq_len(R_REPS), function(r)
    one_rep(d$k_studies, d$rho_xy_mean, d$rel_pattern)))
  tr <- d$rho_xy_mean
  data.frame(
    k_studies = d$k_studies, rho_xy_mean = tr, rel_pattern = d$rel_pattern,
    n_reps = nrow(reps),
    coverage_wald = mean(reps[, "wald_lo"] <= tr & tr <= reps[, "wald_hi"]),
    coverage_knha = mean(reps[, "knha_lo"] <= tr & tr <= reps[, "knha_hi"]),
    width_wald = mean(reps[, "wald_hi"] - reps[, "wald_lo"]),
    width_knha = mean(reps[, "knha_hi"] - reps[, "knha_lo"])
  )
})
res <- do.call(rbind, out)

dir.create("analysis/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "analysis/results/substudy_knapp_hartung.rds")
write.csv(res, "analysis/results/substudy_knapp_hartung.csv", row.names = FALSE)
print(res, row.names = FALSE)
