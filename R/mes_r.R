# mes_r.R: correlation corrected for measurement error, with Fisher-z CI


#' Correlation corrected for measurement error
#'
#' @description
#' Pearson correlation disattenuated for unreliability in one or both
#' variables, with a Fisher-z confidence interval. The interval is built
#' on the observed correlation and its endpoints are then divided by the
#' same attenuation factor as the point estimate, which is exact when the
#' reliabilities are treated as known constants.
#'
#' @details
#' Point correction (Spearman, 1904):
#' \deqn{r^* = \frac{r_{xy}}{\sqrt{r_{xx} \, r_{yy}}}}
#' When \code{ryy} is omitted, only one variable is treated as fallible
#' and the divisor is \eqn{\sqrt{r_{xx}}} (a warning makes the assumption
#' explicit).
#'
#' Interval: with \eqn{z = \mathrm{atanh}(r_{xy})} and standard error
#' \eqn{1/\sqrt{n - 3}}, the Fisher-z interval for the observed
#' correlation is back-transformed and its endpoints divided by the
#' attenuation factor. Because the factor is a fixed positive constant
#' under the known-reliability assumption, the corrected interval has the
#' same coverage as the observed one. Endpoints are truncated to
#' \eqn{[-1, 1]} when the correction pushes them outside the admissible
#' range.
#'
#' @param r Observed Pearson correlation (single value in \eqn{(-1, 1)}).
#' @param n Sample size the correlation is based on (must be at least 4).
#' @param rxx Reliability of the first variable (\eqn{r_{xx} \in (0, 1]}).
#' @param ryy Reliability of the second variable. Omit for the one-sided
#'   correction.
#' @param ci Logical. If TRUE (default), the Fisher-z interval is
#'   computed.
#' @param conf_level Confidence level (default: 0.95).
#'
#' @return An object of class 'mes' with \code{type = "r"}.
#'
#' @note
#' The Fisher-z interval treats \code{rxx} and \code{ryy} as known
#' constants. Unlike \code{\link{mes_d}} and \code{\link{mes_g}}, this
#' function does not (in this version) accept a standard error for the
#' reliabilities, so it does not propagate the extra uncertainty that
#' arises when a reliability is itself a sample estimate. Read the
#' interval as approximate whenever \code{rxx} or \code{ryy} comes from a
#' sample rather than a known population value.
#'
#' @references
#' Spearman, C. (1904). The proof and measurement of association between
#' two things. \emph{American Journal of Psychology}, 15, 72--101.
#'
#' Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis:
#' Correcting error and bias in research findings} (2nd ed.). Sage.
#'
#' @examples
#' # Both variables fallible
#' mes_r(r = 0.30, n = 120, rxx = 0.85, ryy = 0.78)
#'
#' # One-sided correction (second variable treated as error-free)
#' mes_r(r = 0.30, n = 120, rxx = 0.85)
#'
#' @export
mes_r <- function(r,
                  n,
                  rxx,
                  ryy = NULL,
                  ci = TRUE,
                  conf_level = 0.95) {

  call_ <- match.call()

  if (length(r) != 1L || !is.numeric(r) || is.na(r) || abs(r) >= 1) {
    stop("'r' must be a single correlation strictly between -1 and 1.",
         call. = FALSE)
  }
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || n < 4) {
    stop("'n' must be a single sample size of at least 4.", call. = FALSE)
  }
  .check_conf_level(conf_level)

  # validation of rxx/ryy (and the one-sided warning) happens here
  r_corrected <- disattenuate(r, rxx = rxx, ryy = ryy, type = "r")
  atten <- if (is.null(ryy)) sqrt(rxx) else sqrt(rxx * ryy)

  ci_vals <- NULL
  if (ci) {
    z    <- atanh(r)
    se_z <- 1 / sqrt(n - 3)
    crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    ci_obs  <- tanh(c(z - crit * se_z, z + crit * se_z))
    ci_vals <- pmax(pmin(ci_obs / atten, 1), -1)
  }

  .new_mes(
    estimate     = r_corrected,
    raw_estimate = r,
    ci           = ci_vals,
    ci_method    = if (ci) "fisher_z" else NA_character_,
    conf_level   = conf_level,
    type         = "r",
    reliability  = if (is.null(ryy)) list(rxx = rxx)
                   else list(rxx = rxx, ryy = ryy),
    n            = list(total = as.integer(n)),
    call         = call_
  )
}
