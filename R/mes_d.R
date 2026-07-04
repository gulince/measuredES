# mes_d.R: Cohen's d with reliability correction


#' Cohen's d corrected for measurement error
#'
#' @description
#' Cohen's d for two independent groups, with Spearman's (1904) correction
#' for the reliability of the outcome. Pass raw vectors, summary statistics,
#' or a \code{y ~ group} formula, whichever you happen to have.
#'
#' @details
#' Classical Cohen's d:
#' \deqn{d = \frac{\bar{X}_1 - \bar{X}_2}{S_p}}
#' where \eqn{S_p} is the pooled standard deviation.
#'
#' Corrected estimate:
#' \deqn{d^* = \frac{d}{\sqrt{r_{xx}}}}
#'
#' \strong{What the corrected value means:} it estimates what the
#' standardized mean difference would look like if your outcome were
#' measured perfectly. That interpretation rests on the parallel-tests
#' assumption (Lord & Novick, 1968); when reliability is itself estimated
#' the resulting interval is narrower than it should be (Charles, 2005),
#' so read the bounds as approximate.
#'
#' @param x1 Numeric vector for group 1 (raw-data path).
#' @param x2 Numeric vector for group 2 (raw-data path).
#' @param rxx Reliability of the dependent variable (\eqn{r_{xx} \in (0, 1]}).
#' @param rxx_se Optional standard error of the reliability estimate. When
#'   supplied, the analytical intervals (Hunter-Schmidt, Bonett) add the
#'   first-order delta-method term
#'   \eqn{d^2 \, \mathrm{Var}(\hat{r}_{xx}) / (4 r_{xx}^3)} to the variance,
#'   so the interval reflects that the reliability itself was estimated.
#'   Ignored (with a warning) for bootstrap intervals, which resample the
#'   cases but hold \code{rxx} fixed.
#' @param means Length-2 numeric vector; \code{c(M1, M2)} (summary path).
#' @param sds Length-2 numeric vector; \code{c(SD1, SD2)} (summary path).
#' @param ns Length-2 integer vector; \code{c(n1, n2)} (summary path).
#' @param formula A two-sided formula \code{outcome ~ group} where \code{group}
#'   has exactly two levels (formula path). Requires \code{data}.
#' @param data A data frame containing the variables in \code{formula}.
#' @param ci Logical. If TRUE, a confidence interval is computed
#'   (default: TRUE).
#' @param ci_method CI method: \code{"hunter_schmidt"} (analytical),
#'   \code{"bootstrap"} (resampling), \code{"bonett"} (Bonett, 2008
#'   refinement). Bootstrap is available only when raw data are supplied.
#' @param conf_level Confidence level (default: 0.95).
#' @param n_boot Number of bootstrap replications (default: 2000).
#'
#' @return An object of class 'mes'. Components include:
#' \describe{
#'   \item{estimate}{Corrected Cohen's d}
#'   \item{raw_estimate}{Uncorrected (observed) Cohen's d}
#'   \item{ci}{Confidence interval}
#'   \item{reliability}{Reliability information used}
#'   \item{n}{Group sizes}
#' }
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#' sciences} (2nd ed.). Lawrence Erlbaum.
#'
#' Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis:
#' Correcting error and bias in research findings} (2nd ed.). Sage.
#'
#' Bonett, D. G. (2008). Confidence intervals for standardized linear
#' contrasts of means. \emph{Psychological Methods}, 13(2), 99--109.
#'
#' Lord, F. M., & Novick, M. R. (1968). \emph{Statistical theories of
#' mental test scores}. Addison-Wesley.
#'
#' @examples
#' # --- Example 1: From raw data ---
#' set.seed(2026)
#' group1 <- rnorm(50, mean = 100, sd = 15)
#' group2 <- rnorm(50, mean = 110, sd = 15)
#'
#' result <- mes_d(x1 = group1, x2 = group2, rxx = 0.85)
#' print(result)
#'
#' # --- Example 2: From summary statistics ---
#' mes_d(means = c(100, 108),
#'       sds = c(15, 14),
#'       ns = c(80, 75),
#'       rxx = 0.78)
#'
#' # --- Example 3: Formula interface ---
#' df <- data.frame(score = c(group1, group2),
#'                  group = rep(c("control", "treatment"), each = 50))
#' mes_d(formula = score ~ group, data = df, rxx = 0.85)
#'
#' # --- Example 4: Bootstrap confidence interval ---
#' \donttest{
#' mes_d(x1 = group1, x2 = group2, rxx = 0.85,
#'       ci_method = "bootstrap", n_boot = 1000)
#' }
#'
#' @export
mes_d <- function(x1 = NULL,
                  x2 = NULL,
                  rxx,
                  rxx_se = NULL,
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

  if (!is.null(rxx_se)) {
    if (length(rxx_se) != 1L || !is.numeric(rxx_se) ||
        is.na(rxx_se) || rxx_se <= 0) {
      stop("'rxx_se' must be a single positive number.", call. = FALSE)
    }
  }

  # strict = FALSE: disattenuate() below already warns on very low rxx
  .check_reliability(rxx, strict = FALSE)

  if (!is.null(formula)) {
    res <- .resolve_formula(formula, data)
    x1  <- res$x1
    x2  <- res$x2
  }

  raw_data_path <- !is.null(x1) && !is.null(x2)
  summary_path  <- !is.null(means) && !is.null(sds) && !is.null(ns)

  if (!raw_data_path && !summary_path) {
    stop("Provide (x1, x2), (means, sds, ns), or (formula, data).",
         call. = FALSE)
  }

  if (raw_data_path && summary_path) {
    warning("Both raw data and summary statistics provided; raw data used.",
            call. = FALSE)
    summary_path <- FALSE
  }

  if (raw_data_path) {

    x1 <- .check_numeric_vector(x1, "x1")
    x2 <- .check_numeric_vector(x2, "x2")

    M1 <- mean(x1); M2 <- mean(x2)
    s1 <- stats::sd(x1); s2 <- stats::sd(x2)
    n1 <- length(x1); n2 <- length(x2)
    Sp <- .pooled_sd(x1, x2)

  } else {
    if (length(means) != 2L || length(sds) != 2L || length(ns) != 2L) {
      stop("'means', 'sds', and 'ns' must be length-2 vectors.",
           call. = FALSE)
    }
    if (!is.numeric(means) || !is.numeric(sds) || !is.numeric(ns)) {
      stop("'means', 'sds', and 'ns' must be numeric.", call. = FALSE)
    }
    if (any(sds <= 0)) stop("Standard deviations must be positive.",
                            call. = FALSE)
    if (any(ns < 2L))  stop("Each group must have at least 2 observations.",
                            call. = FALSE)

    M1 <- means[1L]; M2 <- means[2L]
    s1 <- sds[1L];   s2 <- sds[2L]
    n1 <- ns[1L];    n2 <- ns[2L]
    Sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  }

  if (Sp == 0) {
    stop("Pooled standard deviation is zero; effect size is undefined.",
         call. = FALSE)
  }

  d_raw <- (M1 - M2) / Sp

  if (!ci) {
    ci_vals <- NULL
  } else if (ci_method == "bootstrap") {
    if (!is.null(rxx_se)) {
      warning("'rxx_se' is ignored for bootstrap intervals: resampling ",
              "holds rxx fixed. Use an analytical ci_method to propagate ",
              "reliability uncertainty.", call. = FALSE)
    }
    if (raw_data_path) {
      ci_vals <- .boot_ci_d(x1, x2, rxx,
                            conf_level = conf_level, n_boot = n_boot)
    } else {
      warning("Bootstrap requires raw data; falling back to Hunter-Schmidt.",
              call. = FALSE)
      ci_method <- "hunter_schmidt"
      ci_vals   <- .hs_ci_d(d_raw, rxx, n1, n2, conf_level, rxx_se)
    }
  } else if (ci_method == "hunter_schmidt") {
    ci_vals <- .hs_ci_d(d_raw, rxx, n1, n2, conf_level, rxx_se)
  } else if (ci_method == "bonett") {
    ci_vals <- .bonett_ci_d(d_raw, rxx, n1, n2, s1, s2, conf_level, rxx_se)
  }

  d_corrected <- disattenuate(d_raw, rxx = rxx, type = "d")

  .new_mes(
    estimate     = d_corrected,
    raw_estimate = d_raw,
    ci           = ci_vals,
    ci_method    = if (ci) ci_method else NA_character_,
    conf_level   = conf_level,
    type         = "d",
    reliability  = if (is.null(rxx_se)) list(rxx = rxx)
                   else list(rxx = rxx, rxx_se = rxx_se),
    n            = list(n1 = as.integer(n1), n2 = as.integer(n2),
                        total = as.integer(n1 + n2)),
    call         = call_,
    extras       = list(M1 = M1, M2 = M2, s1 = s1, s2 = s2, pooled_sd = Sp)
  )
}


# --- Internal CI helpers ---

#' Hunter-Schmidt analytical CI for corrected d
#'
#' @description
#' Sampling variance from Hedges & Olkin (1985), then divided by rxx
#' per Hunter & Schmidt (2004):
#' \deqn{Var(d) = \frac{n_1 + n_2}{n_1 n_2} + \frac{d^2}{2(n_1 + n_2)}}
#' Treating reliability as fixed, the corrected variance is:
#' \deqn{Var(d^*) = \frac{Var(d)}{r_{xx}}}
#'
#' @section Note:
#' By default reliability is treated as a known constant. If \code{rxx_se}
#' is supplied, the first-order delta-method term
#' \eqn{d^2 \, \mathrm{Var}(\hat{r}_{xx}) / (4 r_{xx}^3)} is added to the
#' variance, which widens the interval to reflect that the reliability was
#' estimated rather than known.
#'
#' @keywords internal
.hs_ci_d <- function(d_raw, rxx, n1, n2, conf_level, rxx_se = NULL) {

  # Sampling variance of observed d (Hedges & Olkin, 1985)
  var_d <- (n1 + n2) / (n1 * n2) + d_raw^2 / (2 * (n1 + n2))

  # Approximate variance of corrected d (Hunter & Schmidt, 2004)
  var_d_corrected <- var_d / rxx
  if (!is.null(rxx_se)) {
    var_d_corrected <- var_d_corrected +
      d_raw^2 * rxx_se^2 / (4 * rxx^3)
  }

  d_corrected <- d_raw / sqrt(rxx)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  se <- sqrt(var_d_corrected)

  c(d_corrected - z * se, d_corrected + z * se)
}


#' Bonett (2008) CI for corrected d
#'
#' @description
#' Bonett's (2008) heteroscedasticity-robust sampling variance for the
#' standardized mean difference. Unlike Hunter-Schmidt / Hedges-Olkin
#' (which assume equal population variances), Bonett's formula uses each
#' group's own variance and yields better coverage when the two groups'
#' SDs differ. The interval uses a t critical value with df = n1+n2-2.
#'
#' Variance (Bonett, 2008, adapted to the corrected d):
#' \deqn{\widehat{Var}(d) =
#'   \frac{d^2 \, (s_1^4/v_1 + s_2^4/v_2)}{2(s_1^2 + s_2^2)^2}
#'   + \frac{s_1^2/v_1 + s_2^2/v_2}{(s_1^2 + s_2^2)/2}}
#' where \eqn{v_i = n_i - 1}. The corrected variance is then \eqn{Var(d)/r_{xx}}.
#'
#' One caveat: Bonett derives this variance for a contrast standardized by
#' \eqn{\sqrt{(s_1^2 + s_2^2)/2}}, whereas \code{mes_d()} standardizes by the
#' pooled SD. The two standardizers are identical when \eqn{n_1 = n_2} and
#' stay close under mild imbalance, so the variance is used here as an
#' approximation for the pooled-SD d.
#'
#' @keywords internal
.bonett_ci_d <- function(d_raw, rxx, n1, n2, s1, s2, conf_level,
                         rxx_se = NULL) {

  df <- n1 + n2 - 2L
  v1 <- n1 - 1L; v2 <- n2 - 1L

  s1sq <- s1^2; s2sq <- s2^2
  sum_sq <- s1sq + s2sq

  term1 <- d_raw^2 * (s1sq^2 / v1 + s2sq^2 / v2) / (2 * sum_sq^2)
  term2 <- (s1sq / v1 + s2sq / v2) / (sum_sq / 2)

  var_d <- term1 + term2
  var_d_corrected <- var_d / rxx
  if (!is.null(rxx_se)) {
    var_d_corrected <- var_d_corrected +
      d_raw^2 * rxx_se^2 / (4 * rxx^3)
  }

  d_corrected <- d_raw / sqrt(rxx)
  t_crit <- stats::qt(1 - (1 - conf_level) / 2, df = df)
  se <- sqrt(var_d_corrected)

  c(d_corrected - t_crit * se, d_corrected + t_crit * se)
}


#' Bootstrap CI for corrected d
#'
#' @description
#' Stratified resampling within each group, computing the corrected d on
#' every replicate. BCa intervals by default; falls back to percentile if
#' BCa fails to converge.
#'
#' @keywords internal
.boot_ci_d <- function(x1, x2, rxx, conf_level, n_boot) {

  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("The 'boot' package is required for bootstrap CIs.",
         call. = FALSE)
  }

  combined <- c(x1, x2)
  group_id <- c(rep(1L, length(x1)), rep(2L, length(x2)))

  stat_fn <- function(data, indices) {
    d_resampled <- data[indices, , drop = FALSE]
    g1 <- d_resampled$value[d_resampled$group == 1L]
    g2 <- d_resampled$value[d_resampled$group == 2L]

    if (length(g1) < 2L || length(g2) < 2L) return(NA_real_)

    Sp <- .pooled_sd(g1, g2)
    if (Sp == 0) return(NA_real_)

    d <- (mean(g1) - mean(g2)) / Sp
    d / sqrt(rxx)
  }

  df_boot <- data.frame(value = combined, group = group_id)
  boot_obj <- boot::boot(df_boot, statistic = stat_fn, R = n_boot,
                         strata = factor(group_id))

  # BCa CI; fall back to percentile, and to NA if both fail
  ci_obj <- tryCatch(
    boot::boot.ci(boot_obj, conf = conf_level, type = "bca"),
    error = function(e) tryCatch(
      boot::boot.ci(boot_obj, conf = conf_level, type = "perc"),
      error = function(e2) NULL))

  if (is.null(ci_obj)) c(NA_real_, NA_real_)
  else if (!is.null(ci_obj$bca)) ci_obj$bca[4:5]
  else if (!is.null(ci_obj$percent)) ci_obj$percent[4:5]
  else c(NA_real_, NA_real_)
}
