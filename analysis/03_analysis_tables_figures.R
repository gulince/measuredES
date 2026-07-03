# 03_analysis_tables_figures.R
# Turns the simulation results into the manuscript tables and figures.

rm(list = ls())

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

results      <- readRDS("analysis/results/simulation_main.rds")
results_meta <- readRDS("analysis/results/simulation_meta.rds")


# ---- Table 1: bias summary (rho_true x delta) ----
table1 <- results |>
  dplyr::filter(reliability_known == TRUE, dist == "normal") |>
  dplyr::group_by(rho_true, delta) |>
  dplyr::summarise(
    bias_obs   = round(mean(bias_obs),  3),
    bias_corr  = round(mean(bias_corr), 3),
    rmse_obs   = round(mean(rmse_obs),  3),
    rmse_corr  = round(mean(rmse_corr), 3),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from  = delta,
    values_from = c(bias_obs, bias_corr, rmse_obs, rmse_corr),
    names_glue  = "{.value}_d{delta}"
  )

write.csv(table1, "analysis/results/table1_bias.csv", row.names = FALSE)
cat("Table 1 saved.\n")


# ---- Table 2: CI coverage ----
table2 <- results |>
  dplyr::filter(dist == "normal") |>
  dplyr::group_by(rho_true, n_per_group, reliability_known) |>
  dplyr::summarise(
    cov_hs    = round(mean(coverage_hs), 3),
    cov_bo    = round(mean(coverage_bo), 3),
    width_hs  = round(mean(width_hs), 3),
    width_bo  = round(mean(width_bo), 3),
    .groups = "drop"
  )

write.csv(table2, "analysis/results/table2_coverage.csv", row.names = FALSE)
cat("Table 2 saved.\n")


# ---- Figure 1: bias (rho_true x delta panels) ----
# facet labels like "rho == 0.60" are rendered as math by label_parsed
fig1_data <- results |>
  dplyr::filter(reliability_known == TRUE, dist == "normal") |>
  dplyr::select(delta, n_per_group, rho_true, bias_obs, bias_corr) |>
  tidyr::pivot_longer(c(bias_obs, bias_corr),
                      names_to = "Estimator", values_to = "Bias") |>
  dplyr::mutate(
    Estimator   = dplyr::recode(Estimator,
                                bias_obs  = "Uncorrected",
                                bias_corr = "Corrected"),
    rho_label   = sprintf("rho == %.2f", rho_true),
    delta_label = sprintf("delta == %.2f", delta)
  )

fig1 <- ggplot(fig1_data,
               aes(x = n_per_group, y = Bias,
                   colour = Estimator, shape = Estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.4) +
  geom_line(aes(group = Estimator), linewidth = 0.4) +
  facet_grid(delta_label ~ rho_label, labeller = label_parsed) +
  scale_x_continuous(breaks = c(30, 60, 120, 240)) +
  scale_colour_manual(values = c("Corrected"   = "black",
                                 "Uncorrected" = "grey55")) +
  labs(x = "Sample size per group", y = "Bias",
       colour = NULL, shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

ggsave("analysis/figures/fig1_bias.png", fig1,
       width = 8, height = 6, dpi = 300)
ggsave("analysis/figures/fig1_bias.pdf", fig1,
       width = 8, height = 6)
cat("Figure 1 saved.\n")


# ---- Figure 2: CI coverage ----
fig2_data <- results |>
  dplyr::filter(dist == "normal") |>
  dplyr::select(delta, n_per_group, rho_true, reliability_known,
                coverage_hs, coverage_bo) |>
  tidyr::pivot_longer(c(coverage_hs, coverage_bo),
                      names_to = "Method", values_to = "Coverage") |>
  dplyr::mutate(
    Method = dplyr::recode(Method,
                           coverage_hs = "Hunter-Schmidt",
                           coverage_bo = "Bonett"),
    Reliability = ifelse(reliability_known,
                         "Reliability known",
                         "Reliability estimated"),
    rho_label   = sprintf("rho == %.2f", rho_true)
  )

fig2 <- ggplot(fig2_data,
               aes(x = n_per_group, y = Coverage,
                   colour = Method, linetype = Reliability)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ rho_label, ncol = 2, labeller = label_parsed) +
  scale_x_continuous(breaks = c(30, 60, 120, 240)) +
  scale_y_continuous(limits = c(0.80, 1.00),
                     breaks = seq(0.80, 1.00, 0.05)) +
  scale_colour_manual(values = c("Hunter-Schmidt" = "grey55",
                                 "Bonett"         = "black")) +
  labs(x = "Sample size per group", y = "Coverage probability",
       colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave("analysis/figures/fig2_coverage.png", fig2,
       width = 8, height = 6, dpi = 300)
ggsave("analysis/figures/fig2_coverage.pdf", fig2,
       width = 8, height = 6)
cat("Figure 2 saved.\n")


# ---- Table 3: meta-analytic comparison ----
table3 <- results_meta |>
  dplyr::select(k_studies, delta_mean, rho_pattern,
                bias_pooled_obs, bias_pooled_corr,
                rmse_pooled_obs, rmse_pooled_corr,
                coverage_obs, coverage_corr) |>
  dplyr::mutate(across(where(is.numeric), \(x) round(x, 3)))

write.csv(table3, "analysis/results/table3_meta.csv", row.names = FALSE)
cat("Table 3 saved.\n")


# ---- Figure 3: meta-analytic bias ----
fig3_data <- results_meta |>
  dplyr::select(k_studies, delta_mean, rho_pattern,
                bias_pooled_obs, bias_pooled_corr) |>
  tidyr::pivot_longer(c(bias_pooled_obs, bias_pooled_corr),
                      names_to = "Estimator", values_to = "Bias") |>
  dplyr::mutate(
    Estimator   = dplyr::recode(Estimator,
                                bias_pooled_obs  = "Uncorrected",
                                bias_pooled_corr = "Corrected"),
    delta_label = sprintf("delta == %.2f", delta_mean),
    rho_label   = ifelse(rho_pattern == "homogeneous",
                         "Homogeneous reliability",
                         "Heterogeneous reliability")
  )

fig3 <- ggplot(fig3_data,
               aes(x = factor(k_studies), y = Bias, fill = Estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  # delta_label is parsed as math, rho_label stays plain text
  facet_grid(delta_label ~ rho_label,
             labeller = labeller(delta_label = label_parsed,
                                 rho_label   = label_value)) +
  scale_fill_manual(values = c("Corrected"   = "grey25",
                               "Uncorrected" = "grey70")) +
  labs(x = "Number of studies (k)", y = "Bias of pooled estimate",
       fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave("analysis/figures/fig3_meta_bias.png", fig3,
       width = 8, height = 5, dpi = 300)
ggsave("analysis/figures/fig3_meta_bias.pdf", fig3,
       width = 8, height = 5)
cat("Figure 3 saved.\n")


cat("\nAll tables and figures generated successfully.\n")
cat("Tables : analysis/results/\n")
cat("Figures: analysis/figures/\n")
