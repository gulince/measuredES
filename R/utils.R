# utils.R: internal helpers: input checks, pooled SD, CI computations


#' Bounds check for a reliability coefficient
#'
#' @description
#' Checks that \eqn{r_{xx}} is a single number in \eqn{(0, 1]}.
#' Zero and negative values make the correction mathematically meaningless,
#' and values above 1 are not valid reliabilities.
#'
#' @param r Reliability coefficient (numeric scalar).
#' @param strict Logical. If TRUE, values close to zero raise a warning
#'   (numerical instability) and \eqn{r = 0} raises an error
#'   (correction is mathematically undefined).
#'
#' @return The validated reliability value (returned invisibly).
#'
#' @keywords internal
.check_reliability <- function(r, strict = TRUE) {

  if (length(r) != 1L) {
    stop("Reliability coefficient must be a single numeric value.",
         call. = FALSE)
  }

  if (!is.numeric(r) || is.na(r)) {
    stop("Reliability coefficient (r) must be numeric and non-missing.",
         call. = FALSE)
  }

  if (r <= 0 || r > 1) {
    stop(sprintf(
      "Reliability must lie in (0, 1]; supplied value: %.4f", r),
      call. = FALSE)
  }

  # Very low reliability makes the correction explode
  if (strict && r < 0.10) {
    warning(sprintf(
      paste0("Reliability is very low (r = %.3f). At this level the ",
             "corrected estimate can be several times the observed value ",
             "and is hard to interpret meaningfully."), r),
      call. = FALSE)
  }

  invisible(r)
}


#' Pooled standard deviation for two independent samples
#'
#' @description
#' Pooled SD using the Hedges & Olkin (1985) formula:
#' \deqn{S_p = \sqrt{\frac{(n_1 - 1)\,s_1^2 + (n_2 - 1)\,s_2^2}{n_1 + n_2 - 2}}}
#'
#' @param x1,x2 Numeric vectors (the two groups).
#'
#' @return The pooled standard deviation (numeric scalar).
#'
#' @keywords internal
.pooled_sd <- function(x1, x2) {

  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]

  n1 <- length(x1)
  n2 <- length(x2)

  if (n1 < 2L || n2 < 2L) {
    stop("Each group must contain at least 2 observations.", call. = FALSE)
  }

  v1 <- stats::var(x1)
  v2 <- stats::var(x2)

  pooled_var <- ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)

  sqrt(pooled_var)
}


#' Hedges' small-sample bias correction factor
#'
#' @description
#' Cohen's d overestimates the population effect in small samples.
#' The Hedges (1981) correction factor is:
#' \deqn{J(\nu) = 1 - \frac{3}{4\nu - 1}}
#' where \eqn{\nu = n_1 + n_2 - 2}.
#'
#' @param df Degrees of freedom.
#'
#' @return The Hedges correction factor J.
#'
#' @keywords internal
.hedges_correction <- function(df) {
  if (df <= 1L) {
    stop("Degrees of freedom must be > 1 for the Hedges correction.",
         call. = FALSE)
  }
  1 - (3 / (4 * df - 1))
}


#' Resolve a y ~ group formula into two numeric vectors
#'
#' @description
#' Shared by \code{mes_d()} and \code{mes_g()} so the formula interface
#' behaves identically in both.
#'
#' @keywords internal
.resolve_formula <- function(formula, data) {
  if (is.null(data)) {
    stop("'data' must be supplied when using 'formula'.", call. = FALSE)
  }
  mf   <- stats::model.frame(formula, data = data)
  if (ncol(mf) != 2L) {
    stop("'formula' must have exactly one grouping variable on the ",
         "right-hand side (outcome ~ group).", call. = FALSE)
  }
  y    <- mf[[1L]]
  grp  <- factor(mf[[2L]])
  lvls <- levels(grp)
  if (length(lvls) != 2L) {
    stop("Formula grouping variable must have exactly 2 levels.",
         call. = FALSE)
  }
  list(x1 = y[grp == lvls[1L]], x2 = y[grp == lvls[2L]])
}


#' Validate a numeric input vector
#'
#' @keywords internal
.check_numeric_vector <- function(x, name = "x", min_n = 2L) {
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be a numeric vector.", name), call. = FALSE)
  }
  x_clean <- x[!is.na(x)]
  if (length(x_clean) < min_n) {
    stop(sprintf("'%s' must contain at least %d complete observations.",
                 name, min_n), call. = FALSE)
  }
  invisible(x_clean)
}


#' Verbal interpretation of effect-size magnitude (Cohen's guidelines)
#'
#' @description
#' Rough Cohen (1988) size label. Vectorized: accepts a numeric vector and
#' returns a character vector of the same length. Treat the cutoffs as
#' ballpark figures, not hard rules: context always matters (Lakens, 2013).
#'
#' @keywords internal
.cohen_interpretation <- function(d, type = "d") {
  if (!type %in% c("d", "g")) return(rep(NA_character_, length(d)))
  ad <- abs(d)
  ifelse(ad < 0.20, "negligible",
    ifelse(ad < 0.50, "small",
      ifelse(ad < 0.80, "medium", "large")))
}


# Argument checks shared by mes_d(), mes_g(), mes_r() and mes_ci().
# Plain comments (not roxygen) on purpose: these helpers need no Rd page.

# Confidence level: a single number strictly between 0 and 1
.check_conf_level <- function(conf_level) {
  if (length(conf_level) != 1L || !is.numeric(conf_level) ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  invisible(conf_level)
}

# Number of bootstrap replications: a single whole number, at least 10
.check_n_boot <- function(n_boot) {
  if (length(n_boot) != 1L || !is.numeric(n_boot) || is.na(n_boot) ||
      n_boot < 10 || n_boot != round(n_boot)) {
    stop("'n_boot' must be a single whole number of at least 10.",
         call. = FALSE)
  }
  invisible(n_boot)
}

# Summary-statistics input: numeric, complete, finite, whole-number ns
.check_summary_stats <- function(means, sds, ns) {
  if (length(means) != 2L || length(sds) != 2L || length(ns) != 2L) {
    stop("'means', 'sds', and 'ns' must be length-2 vectors.",
         call. = FALSE)
  }
  if (!is.numeric(means) || !is.numeric(sds) || !is.numeric(ns)) {
    stop("'means', 'sds', and 'ns' must be numeric.", call. = FALSE)
  }
  if (!all(is.finite(means)) || !all(is.finite(sds)) ||
      !all(is.finite(ns))) {
    stop("'means', 'sds', and 'ns' must be finite and non-missing.",
         call. = FALSE)
  }
  if (any(sds <= 0)) {
    stop("Standard deviations must be positive.", call. = FALSE)
  }
  if (any(ns != round(ns))) {
    stop("'ns' must be whole numbers.", call. = FALSE)
  }
  if (any(ns < 2L)) {
    stop("Each group must have at least 2 observations.", call. = FALSE)
  }
  invisible(NULL)
}


# Suppress R CMD check NOTEs for ggplot2 aes() variable names used in plot.mes()
utils::globalVariables(c("estimate", "label", "lower", "upper", ".data"))
