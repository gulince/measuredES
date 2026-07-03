# 04_example_reanalysis.R
# Worked example for the manuscript (Section 3.4): reanalysis of the
# d-value example dataset distributed with psychmeta. Compares the pooled
# random-effects estimate before and after study-level disattenuation.

suppressMessages({
  library(metafor)
  library(measuredES)
})

data("data_d_meas_multi", package = "psychmeta")
dat <- subset(data_d_meas_multi, construct == "Y")

# Reported d and group sizes reproduce each study exactly through the
# summary-statistic interface (pooled SD = 1, mean difference = d).
studies <- lapply(seq_len(nrow(dat)), function(i) {
  mes_d(means = c(dat$d[i], 0), sds = c(1, 1),
        ns    = c(dat$n1[i], dat$n2[i]),
        rxx   = dat$ryyi[i],
        ci_method = "hunter_schmidt")
})
meta_corr <- mes_meta(studies)

yi_obs <- dat$d
vi_obs <- (dat$n1 + dat$n2) / (dat$n1 * dat$n2) +
          dat$d^2 / (2 * (dat$n1 + dat$n2))

fit_obs  <- rma(yi = yi_obs,       vi = vi_obs,       method = "REML")
fit_corr <- rma(yi = meta_corr$yi, vi = meta_corr$vi, method = "REML")

cat(sprintf("k = %d studies | mean ryy = %.3f (range %.2f-%.2f)\n",
            nrow(dat), mean(dat$ryyi), min(dat$ryyi), max(dat$ryyi)))
cat(sprintf("Uncorrected: pooled d = %.3f, 95%% CI [%.3f, %.3f], tau2 = %.4f\n",
            fit_obs$b, fit_obs$ci.lb, fit_obs$ci.ub, fit_obs$tau2))
cat(sprintf("Corrected  : pooled d = %.3f, 95%% CI [%.3f, %.3f], tau2 = %.4f\n",
            fit_corr$b, fit_corr$ci.lb, fit_corr$ci.ub, fit_corr$tau2))
cat(sprintf("Relative increase in pooled estimate: %.1f%%\n",
            100 * (fit_corr$b / fit_obs$b - 1)))
