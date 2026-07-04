#' measuredES: Measurement-Corrected Effect Sizes
#'
#' @description
#' \pkg{measuredES} computes effect sizes corrected for measurement
#' unreliability. The focus is on the two-group situation: Cohen's \emph{d}
#' and Hedges' \emph{g} under classical test theory (CTT).
#'
#' @section Core functions:
#' \describe{
#'   \item{\code{\link{mes_d}}}{Corrected Cohen's d (raw data, summary stats, or formula)}
#'   \item{\code{\link{mes_g}}}{Corrected Hedges' g (small-sample bias removed first)}
#'   \item{\code{\link{mes_r}}}{Corrected Pearson r with Fisher-z confidence interval}
#'   \item{\code{\link{disattenuate}}}{The bare Spearman correction; accepts vectors}
#'   \item{\code{\link{mes_ci}}}{Swap the CI method on an existing result}
#'   \item{\code{\link{mes_meta}}}{Export to a \code{metafor}-ready data frame}
#'   \item{\code{\link{mes_report}}}{APA-7-style inline text or table}
#' }
#'
#' @section Theoretical foundation:
#' In classical test theory, an observed score \eqn{X = T + E} is the sum of
#' a true score (\eqn{T}) and random error (\eqn{E}). The presence of
#' measurement error attenuates observed standardized effect sizes by a
#' factor of \eqn{\sqrt{r_{xx}}} relative to the true population effect
#' (Spearman, 1904). Hunter and Schmidt (2004) document the implications of
#' this attenuation in meta-analysis in detail.
#'
#' @section Caveats:
#' The correction handles random measurement error only; systematic bias
#' is a different problem entirely. At low reliability (\eqn{r_{xx} < 0.50})
#' the corrected estimate becomes unstable and hard to interpret
#' (Charles, 2005), so treat those results with caution.
#'
#' @author
#' Maintainer: Fazilet Gül İnce Aracı \email{gulincegazi@gmail.com}
#'
#' Authors: Fazilet Gül İnce Aracı, Eren Halil Özberk
#'
#' @references
#' Spearman, C. (1904). The proof and measurement of association between
#' two things. \emph{American Journal of Psychology}, 15, 72--101.
#'
#' Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis:
#' Correcting error and bias in research findings} (2nd ed.). Sage.
#'
#' Charles, E. P. (2005). The correction for attenuation due to measurement
#' error: Clarifying concepts and creating confidence sets.
#' \emph{Psychological Methods}, 10(2), 206--226.
#'
#' Baugh, F. (2002). Correcting effect sizes for score reliability.
#' \emph{Educational and Psychological Measurement}, 62, 254--263.
#'
#' @keywords internal
"_PACKAGE"
