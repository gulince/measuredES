# disattenuate.R: Spearman's (1904) correction for attenuation


#' Correct an effect size for measurement error
#'
#' @description
#' The calculation behind every \code{mes_*()} function. Under classical test
#' theory the observed score is true score plus random error
#' (\eqn{X = T + E}), and that error pulls observed standardized effect
#' sizes toward zero by a factor of \eqn{\sqrt{r_{xx}}}. This function
#' undoes that shrinkage given a reliability estimate.
#'
#' @details
#' For standardized mean differences (Cohen's d, Hedges' g), the correction
#' is:
#' \deqn{d_{true} = \frac{d_{observed}}{\sqrt{r_{xx}}}}
#'
#' For bivariate correlations (Pearson r, partial r), both variables can be
#' unreliable:
#' \deqn{r_{true} = \frac{r_{xy}}{\sqrt{r_{xx} \cdot r_{yy}}}}
#' When \code{ryy} is omitted, only the dependent variable (or criterion) is
#' treated as unreliable; the function issues a warning in that case so the
#' assumption is explicit.
#'
#' \strong{Estimated reliabilities (Charles, 2005; Hunter & Schmidt, 2004):}
#' when the reliability is itself a sample estimate, the sampling
#' distribution of the corrected effect is wider than an interval that
#' treats the reliability as fixed. \code{mes_d()} and \code{mes_g()} add a
#' first-order delta-method term when a standard error is supplied through
#' \code{rxx_se}; \code{mes_r()} and the bootstrap intervals treat the
#' reliabilities as known constants.
#'
#' @param effect Observed effect size. May be a single value or a numeric
#'   vector (e.g., a column of \emph{d} values from a study database).
#'   \code{NA} elements pass through as \code{NA}.
#' @param rxx Reliability of the dependent variable (or first variable when
#'   \code{type = "r"}). Must be a single value in \eqn{(0, 1]}.
#' @param ryy Reliability of the second variable; used only when
#'   \code{type = "r"} for a two-sided correction. Omit for the one-sided
#'   correction.
#' @param type Effect-size type:
#'   \itemize{
#'     \item \code{"d"}: Standardized mean difference (Cohen's d / Hedges' g).
#'       Divides by \eqn{\sqrt{r_{xx}}}.
#'     \item \code{"r"}: Pearson correlation. Divides by
#'       \eqn{\sqrt{r_{xx} \cdot r_{yy}}}. For a full workflow that also
#'       returns a Fisher-z confidence interval, use \code{\link{mes_r}}.
#'   }
#'
#' @return A numeric vector of corrected effect sizes with the same length
#'   as \code{effect}.
#'
#' @references
#' Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis:
#' Correcting error and bias in research findings} (2nd ed.). Sage.
#'
#' Spearman, C. (1904). The proof and measurement of association between two
#' things. \emph{American Journal of Psychology}, 15, 72--101.
#'
#' Charles, E. P. (2005). The correction for attenuation due to measurement
#' error: Clarifying concepts and creating confidence sets.
#' \emph{Psychological Methods}, 10(2), 206--226.
#'
#' Baugh, F. (2002). Correcting effect sizes for score reliability: A reminder
#' that measurement and substantive issues are linked inextricably.
#' \emph{Educational and Psychological Measurement}, 62, 254--263.
#'
#' @examples
#' # Single value: Cohen's d at rxx = 0.80
#' disattenuate(effect = 0.50, rxx = 0.80, type = "d")
#'
#' # Vector of d values (e.g., from a study database column)
#' disattenuate(effect = c(0.20, 0.45, 0.63, 0.80), rxx = 0.75, type = "d")
#'
#' # Two-sided correction for a Pearson r
#' disattenuate(effect = 0.30, rxx = 0.85, ryy = 0.78, type = "r")
#'
#' @export
disattenuate <- function(effect,
                         rxx,
                         ryy = NULL,
                         type = c("d", "r")) {

  type <- match.arg(type)

  if (!is.numeric(effect)) {
    stop("'effect' must be a numeric vector.", call. = FALSE)
  }

  .check_reliability(rxx)

  if (type == "d") {
    if (!is.null(ryy)) {
      warning("'ryy' is ignored when type = \"d\"; the d-family ",
              "correction uses only the outcome reliability 'rxx'.",
              call. = FALSE)
    }
    effect / sqrt(rxx)

  } else {
    if (is.null(ryy)) {
      warning("'ryy' not supplied; applying a one-sided correction that ",
              "treats the second variable as measured without error. ",
              "If both variables are fallible, pass 'ryy' explicitly.",
              call. = FALSE)
      effect / sqrt(rxx)
    } else {
      .check_reliability(ryy)
      corrected_r <- effect / sqrt(rxx * ryy)
      out_of_bounds <- !is.na(corrected_r) & abs(corrected_r) > 1
      if (any(out_of_bounds)) {
        warning(sprintf(
          "%d corrected correlation(s) have |r*| > 1, likely caused by sampling error in the reliability estimates or model misspecification (Hunter & Schmidt, 2004).",
          sum(out_of_bounds)), call. = FALSE)
      }
      corrected_r
    }
  }
}
