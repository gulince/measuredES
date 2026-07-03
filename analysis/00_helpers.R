# 00_helpers.R: data generation and effect-size computation for the simulation

# Parameterization (X = T + E):
#   Var(T) = 1, Var(E) = (1 - rho)/rho, so Var(X) = 1/rho and
#   Var(T)/Var(X) = rho. delta is defined on T; on X it attenuates to
#   d_obs = delta * sqrt(rho), hence d_corrected = d_obs / sqrt(rho).

generate_two_groups <- function(n1, n2, delta, rho, dist = c("normal", "skewed")) {

  dist <- match.arg(dist)

  sigma_E <- sqrt((1 - rho) / rho)

  if (dist == "normal") {
    T1 <- rnorm(n1, mean = 0,     sd = 1)
    T2 <- rnorm(n2, mean = delta, sd = 1)
    E1 <- rnorm(n1, mean = 0, sd = sigma_E)
    E2 <- rnorm(n2, mean = 0, sd = sigma_E)
  } else if (dist == "skewed") {
    # standardized chi-sq(4): skewness = sqrt(8/4) ~ 1.41
    skew_sd_factor <- sqrt(2 * 4)
    T1 <- (rchisq(n1, df = 4) - 4) / skew_sd_factor
    T2 <- (rchisq(n2, df = 4) - 4) / skew_sd_factor + delta
    E1 <- (rchisq(n1, df = 4) - 4) / skew_sd_factor * sigma_E
    E2 <- (rchisq(n2, df = 4) - 4) / skew_sd_factor * sigma_E
  }

  X1 <- T1 + E1
  X2 <- T2 + E2

  list(X1 = X1, X2 = X2, T1 = T1, T2 = T2)
}


# Item-level data for the estimated-reliability conditions: a k-item
# tau-equivalent test (X_ij = T_i + e_ij) so Cronbach's alpha can be computed.
generate_test_data <- function(n1, n2, delta, rho, n_items = 10, dist = "normal") {

  # per-item reliability chosen so the total test hits rho (Spearman-Brown,
  # solved for rho_item)
  rho_item <- rho / (k_items_to_factor(n_items, rho))

  if (dist == "normal") {
    T1 <- rnorm(n1, 0, 1)
    T2 <- rnorm(n2, delta, 1)
  } else {
    sd_factor <- sqrt(8)
    T1 <- (rchisq(n1, 4) - 4) / sd_factor
    T2 <- (rchisq(n2, 4) - 4) / sd_factor + delta
  }

  T_all <- c(T1, T2)
  N <- n1 + n2

  sigma_e_item <- sqrt((1 - rho_item) / rho_item)

  E_mat <- matrix(rnorm(N * n_items, 0, sigma_e_item),
                  nrow = N, ncol = n_items)
  X_mat <- T_all + E_mat   # each item = T + e_ij

  X_total <- rowSums(X_mat)

  list(
    X_mat   = X_mat,
    X_total = X_total,
    group   = c(rep(1L, n1), rep(2L, n2)),
    T_all   = T_all
  )
}

# Spearman-Brown denominator for the per-item rho
k_items_to_factor <- function(k, rho_target) {
  k - (k - 1) * rho_target
}


cronbach_alpha <- function(item_matrix) {
  k <- ncol(item_matrix)
  item_vars <- apply(item_matrix, 2, stats::var)
  total_var <- stats::var(rowSums(item_matrix))
  alpha <- (k / (k - 1)) * (1 - sum(item_vars) / total_var)
  alpha
}


pooled_sd_fast <- function(x1, x2) {
  n1 <- length(x1); n2 <- length(x2)
  v1 <- stats::var(x1); v2 <- stats::var(x2)
  sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
}


one_replication <- function(n1, n2, delta, rho_true,
                            reliability_known = TRUE,
                            dist = "normal",
                            n_items = 10) {

  if (reliability_known) {
    dat <- generate_two_groups(n1, n2, delta, rho_true, dist)
    X1 <- dat$X1; X2 <- dat$X2
    rho_hat <- rho_true
  } else {
    dat <- generate_test_data(n1, n2, delta, rho_true, n_items, dist)
    X_total <- dat$X_total
    X1 <- X_total[dat$group == 1L]
    X2 <- X_total[dat$group == 2L]
    rho_hat <- cronbach_alpha(dat$X_mat)
    rho_hat <- max(min(rho_hat, 0.999), 0.05)  # sample alpha can be <= 0
  }

  Sp <- pooled_sd_fast(X1, X2)
  if (Sp == 0) return(NULL)

  d_obs       <- (mean(X2) - mean(X1)) / Sp
  d_corrected <- d_obs / sqrt(rho_hat)

  var_d_obs  <- (n1 + n2) / (n1 * n2) + d_obs^2 / (2 * (n1 + n2))
  var_d_corr <- var_d_obs / rho_hat

  # Hunter-Schmidt (z-based)
  z_crit <- stats::qnorm(0.975)
  ci_hs <- d_corrected + c(-1, 1) * z_crit * sqrt(var_d_corr)

  # Bonett (2008) heteroscedasticity-robust variance
  s1sq   <- stats::var(X1); s2sq <- stats::var(X2)
  v1     <- n1 - 1L;        v2   <- n2 - 1L
  sum_sq <- s1sq + s2sq
  term1  <- d_obs^2 * (s1sq^2 / v1 + s2sq^2 / v2) / (2 * sum_sq^2)
  term2  <- (s1sq / v1 + s2sq / v2) / (sum_sq / 2)
  var_d_bonett_corr <- (term1 + term2) / rho_hat
  t_crit <- stats::qt(0.975, df = n1 + n2 - 2L)
  ci_bo  <- d_corrected + c(-1, 1) * t_crit * sqrt(var_d_bonett_corr)

  # bootstrap CI skipped: too expensive at R = 5000 per condition

  data.frame(
    d_obs        = d_obs,
    d_corrected  = d_corrected,
    rho_hat      = rho_hat,
    ci_hs_low    = ci_hs[1L],
    ci_hs_high   = ci_hs[2L],
    ci_bo_low    = ci_bo[1L],
    ci_bo_high   = ci_bo[2L],
    var_d_corr   = var_d_corr
  )
}
