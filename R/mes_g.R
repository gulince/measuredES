# mes_g.R: Hedges' g (small-sample correction first, then reliability)


#' Hedges' g corrected for measurement error
#'
#' @description
#' Hedges' g is just Cohen's d shrunk slightly to remove the small-sample
#' bias, then disattenuated for measurement error. Two steps:
#'
#' **Step 1: Hedges (1981) bias correction.**
#' \deqn{g = J(\nu) \cdot d}
#' where \eqn{J(\nu) = 1 - \frac{3}{4\nu - 1}} and \eqn{\nu = n_1 + n_2 - 2}.
#' J is just under 1 for small samples and approaches 1 as n grows.
#'
#' **Step 2: reliability correction (Hunter & Schmidt, 2004).**
#' \deqn{g^* = \frac{g}{\sqrt{r_{xx}}}}
#'
#' @details
#' Order matters for the variance. Since g = J·d, the sampling variance
#' picks up a J^2 factor before the reliability correction is applied:
#' \deqn{\text{Var}(g^*) = \frac{J^2 \cdot \text{Var}(d)}{r_{xx}}}
#' (Hunter & Schmidt, 2004, p. 99; Bonett, 2008).
#'
#' @inheritParams mes_d
#'
#' @return An object of class 'mes'.
#'
#' @references
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#' size and related estimators. \emph{Journal of Educational Statistics},
#' 6(2), 107--128.
#'
#' Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for
#' meta-analysis}. Academic Press.
#'
#' Lakens, D. (2013). Calculating and reporting effect sizes to facilitate
#' cumulative science: A practical primer for t-tests and ANOVAs.
#' \emph{Frontiers in Psychology}, 4, 863.
#'
#' @examples
#' set.seed(2026)
#' control <- rnorm(25, 50, 10)
#' treat   <- rnorm(25, 56, 10)
#'
#' # Small sample: Hedges' g is less biased than Cohen's d here
#' mes_g(x1 = control, x2 = treat, rxx = 0.82)
#'
#' # Same call via the formula interface
#' df <- data.frame(score = c(control, treat),
#'                  group = rep(c("control", "treat"), each = 25))
#' mes_g(formula = score ~ group, data = df, rxx = 0.82)
#'
#' @export
mes_g <- function(x1 = NULL,
                  x2 = NULL,
                  rxx,
                  means = NULL,
                  sds = NULL,
                  ns = NULL,
                  formula = NULL,
                  data = NULL,
                  ci = TRUE,
                  ci_method = c("hunter_schmidt", "bootstrap", "bonett"),
                  conf_level = 0.95,
                  n_boot = 2000L) {

  ci_method <- match.arg(ci_method)
  call_ <- match.call()

  if (conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be in (0, 1).", call. = FALSE)
  }

  # resolve the formula early: the bootstrap path needs x1/x2
  if (!is.null(formula)) {
    res <- .resolve_formula(formula, data)
    x1  <- res$x1
    x2  <- res$x2
  }

  # the bootstrap resamples x1/x2 directly, so drop NAs up front
  if (!is.null(x1) && !is.null(x2)) {
    x1 <- .check_numeric_vector(x1, "x1")
    x2 <- .check_numeric_vector(x2, "x2")
  }

  d_obj <- mes_d(x1 = x1, x2 = x2, rxx = rxx,
                 means = means, sds = sds, ns = ns,
                 ci = FALSE,
                 conf_level = conf_level)

  n1 <- d_obj$n$n1
  n2 <- d_obj$n$n2
  df <- n1 + n2 - 2L
  s1 <- d_obj$extras$s1
  s2 <- d_obj$extras$s2

  J  <- .hedges_correction(df)

  g_raw       <- J * d_obj$raw_estimate
  g_corrected <- g_raw / sqrt(rxx)

  # Var(g*) = J^2 * Var(d) / rxx
  ci_vals <- NULL
  if (ci) {
    if (ci_method == "bootstrap" && !is.null(x1) && !is.null(x2)) {
      ci_vals <- .boot_ci_g(x1, x2, rxx, J, conf_level, n_boot)
    } else {
      if (ci_method == "bootstrap") {
        warning("Bootstrap requires raw data; falling back to Hunter-Schmidt.",
                call. = FALSE)
        ci_method <- "hunter_schmidt"
      }

      if (ci_method == "bonett") {
        v1 <- n1 - 1L; v2 <- n2 - 1L
        s1sq <- s1^2; s2sq <- s2^2; sum_sq <- s1sq + s2sq
        var_d <- d_obj$raw_estimate^2 *
                   (s1sq^2 / v1 + s2sq^2 / v2) / (2 * sum_sq^2) +
                 (s1sq / v1 + s2sq / v2) / (sum_sq / 2)
        crit <- stats::qt(1 - (1 - conf_level) / 2, df = df)
      } else {
        # Hunter-Schmidt / Hedges & Olkin asymptotic variance
        var_d <- (n1 + n2) / (n1 * n2) +
                 d_obj$raw_estimate^2 / (2 * (n1 + n2))
        crit <- stats::qnorm(1 - (1 - conf_level) / 2)
      }

      var_g_corrected <- (J^2) * var_d / rxx
      se <- sqrt(var_g_corrected)
      ci_vals <- c(g_corrected - crit * se, g_corrected + crit * se)
    }
  }

  .new_mes(
    estimate     = g_corrected,
    raw_estimate = g_raw,
    ci           = ci_vals,
    ci_method    = if (ci) ci_method else NA_character_,
    conf_level   = conf_level,
    type         = "g",
    reliability  = list(rxx = rxx),
    n            = list(n1 = as.integer(n1), n2 = as.integer(n2),
                        total = as.integer(n1 + n2)),
    call         = call_,
    extras       = list(hedges_J    = J,
                        cohen_d_raw = d_obj$raw_estimate,
                        s1 = s1, s2 = s2)
  )
}


#' Bootstrap CI for g
#'
#' @keywords internal
.boot_ci_g <- function(x1, x2, rxx, J, conf_level, n_boot) {

  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("The 'boot' package is required for bootstrap CIs.",
         call. = FALSE)
  }

  combined <- c(x1, x2)
  group_id <- c(rep(1L, length(x1)), rep(2L, length(x2)))

  stat_fn <- function(data, indices) {
    rs <- data[indices, , drop = FALSE]
    g1 <- rs$value[rs$group == 1L]
    g2 <- rs$value[rs$group == 2L]
    if (length(g1) < 2L || length(g2) < 2L) return(NA_real_)
    Sp <- .pooled_sd(g1, g2)
    if (Sp == 0) return(NA_real_)
    d <- (mean(g1) - mean(g2)) / Sp
    (J * d) / sqrt(rxx)
  }

  df_boot <- data.frame(value = combined, group = group_id)
  boot_obj <- boot::boot(df_boot, statistic = stat_fn, R = n_boot,
                         strata = factor(group_id))

  ci_obj <- tryCatch(
    boot::boot.ci(boot_obj, conf = conf_level, type = "bca"),
    error = function(e) tryCatch(
      boot::boot.ci(boot_obj, conf = conf_level, type = "perc"),
      error = function(e2) NULL))

  if (is.null(ci_obj)) c(NA_real_, NA_real_)
  else if (!is.null(ci_obj$bca))     ci_obj$bca[4:5]
  else if (!is.null(ci_obj$percent)) ci_obj$percent[4:5]
  else c(NA_real_, NA_real_)
}
