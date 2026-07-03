# disattenuate.R: Spearman's (1904) correction for attenuation


#' Correct an effect size for measurement error
#'
#' @description
#' The math behind every \code{mes_*()} function. Under classical test
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
#' \strong{One thing to keep in mind (Charles, 2005; Hunter & Schmidt, 2004):}
#' if your reliability is itself a sample estimate (most of the time it
#' is), the corrected effect's sampling distribution is wider than the
#' analytical formulas in this package suggest. All three CI methods here
#' treat \code{rxx} as a known constant; for the full uncertainty you'd
#' need a Bayesian model or a bootstrap that resamples reliability too.
#'
#' @param effect Observed effect-size value (numeric scalar).
#' @param rxx Reliability of the dependent variable (or of the single scale),
#'   \eqn{r_{xx} \in (0, 1]}.
#' @param ryy Reliability of the second variable (used only when
#'   \code{type = "r"} for two-variable correction).
#' @param type Effect-size type:
#'   \itemize{
#'     \item \code{"d"}: Standardized mean difference (Cohen's d / Hedges' g).
#'       Divides by \eqn{\sqrt{r_{xx}}}.
#'     \item \code{"r"}: Correlation coefficient. Divides by
#'       \eqn{\sqrt{r_{xx} \cdot r_{yy}}}. \code{measuredES} doesn't ship a
#'       dedicated correlation wrapper; for a full workflow see
#'       \code{psych::correct.cor()}.
#'   }
#'
#' @return The corrected effect size (numeric scalar).
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
#' # Reliability correction for Cohen's d (alpha = 0.80)
#' disattenuate(effect = 0.50, rxx = 0.80, type = "d")
#'
#' # Two-sided correction for Pearson r
#' disattenuate(effect = 0.30, rxx = 0.85, ryy = 0.78, type = "r")
#'
#' @export
disattenuate <- function(effect,
                         rxx,
                         ryy = NULL,
                         type = c("d", "r")) {

  type <- match.arg(type)

  if (length(effect) != 1L || !is.numeric(effect) || is.na(effect)) {
    stop("'effect' must be a single numeric value.", call. = FALSE)
  }

  .check_reliability(rxx)

  corrected <- switch(
    type,

    "d" = {
      if (!is.null(ryy)) {
        warning("'ryy' is ignored when type = \"d\"; the d-family ",
                "correction uses only the outcome reliability 'rxx'.",
                call. = FALSE)
      }
      effect / sqrt(rxx)
    },

    "r" = {
      if (is.null(ryy)) {
        warning("'ryy' not supplied; applying a one-sided correction that ",
                "treats the second variable as measured without error. ",
                "If both variables are fallible, pass 'ryy' explicitly.",
                call. = FALSE)
        effect / sqrt(rxx)
      } else {
        .check_reliability(ryy)
        corrected_r <- effect / sqrt(rxx * ryy)
        if (abs(corrected_r) > 1) {
          warning(sprintf(
            paste0("Corrected correlation |r| > 1 (r = %.3f). Likely caused ",
                   "by sampling error in the reliability estimates or model ",
                   "misspecification (Hunter & Schmidt, 2004, p. 96)."),
            corrected_r), call. = FALSE)
        }
        corrected_r
      }
    }
  )

  corrected
}


#' Quick check of a reliability coefficient
#'
#' @description
#' Prints a quick verdict on a reliability value (excellent, good,
#' acceptable, etc.) along with how much score variance is error and
#' how much an observed effect would be shrunk by it. Handy as a
#' sanity check before running the analysis.
#'
#' @param r Reliability coefficient.
#'
#' @return Invisibly, the supplied \code{r}; prints messages.
#'
#' @examples
#' reliability_check(0.85)
#'
#' @export
reliability_check <- function(r) {
  .check_reliability(r, strict = TRUE)

  level <- if (r >= 0.90)      "excellent"
           else if (r >= 0.80) "good"
           else if (r >= 0.70) "acceptable"
           else if (r >= 0.60) "questionable"
           else                "inadequate"

  message(sprintf("r_xx = %.3f  (%s)", r, level))
  message(sprintf("About %.1f%% of score variance is error.", (1 - r) * 100))
  message(sprintf("Observed effects are about %.1f%% of the true effect (sqrt(r) = %.3f).",
                  sqrt(r) * 100, sqrt(r)))

  invisible(r)
}
