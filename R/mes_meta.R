# mes_meta.R: turn mes objects into a data frame metafor can use


#' Convert 'mes' objects to a meta-analysis-ready format
#'
#' @description
#' Transforms one or more 'mes' objects into a data frame that can be passed
#' directly to \code{metafor::rma()} or used in Hunter--Schmidt psychometric
#' meta-analysis procedures.
#'
#' @details
#' The sampling variance \code{vi} is returned on the same scale as the
#' corrected effect \code{yi}. For \emph{d}- and \emph{g}-type objects it
#' comes from the Hedges--Olkin (1985) large-sample formula divided by the
#' reliability (with the delta-method term added when \code{rxx_se} was
#' supplied). For correlations it is
#' \deqn{\widehat{Var}(r^*) =
#'   \frac{(1 - r^2)^2}{(n - 1)\, r_{xx} r_{yy}},}
#' the large-sample variance of the observed correlation scaled by the
#' (fixed) attenuation factor; one-sided corrections drop \eqn{r_{yy}}.
#' Symmetric analytical intervals (Hunter--Schmidt, Bonett) are inverted
#' back to a variance directly, whereas asymmetric BCa bootstrap intervals
#' fall back to the analytical formula.
#'
#' @param ... One or more 'mes' objects (or a list of 'mes' objects).
#' @param study_id Optional study identifiers (character vector).
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{study}{Study identifier}
#'   \item{yi}{Corrected effect size}
#'   \item{vi}{Approximate sampling variance of the effect size}
#'   \item{sei}{Standard error}
#'   \item{ni}{Total sample size}
#'   \item{rxx}{Reliability of the dependent variable (or of the first
#'     variable for correlations)}
#'   \item{type}{Effect-size type}
#' }
#'
#' @references
#' Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis}
#' (2nd ed.). Sage.
#'
#' Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor
#' package. \emph{Journal of Statistical Software}, 36(3), 1--48.
#'
#' @examples
#' set.seed(2026)
#' c1 <- rnorm(40, 100, 15); t1 <- rnorm(40, 108, 15)
#' c2 <- rnorm(60, 100, 15); t2 <- rnorm(60, 105, 15)
#'
#' s1 <- mes_d(c1, t1, rxx = 0.85)
#' s2 <- mes_d(c2, t2, rxx = 0.78)
#'
#' meta_df <- mes_meta(s1, s2, study_id = c("Study_2024", "Study_2025"))
#' meta_df
#'
#' # Corrected correlations combine the same way
#' r1 <- mes_r(r = 0.30, n = 120, rxx = 0.85, ryy = 0.78)
#' r2 <- mes_r(r = 0.42, n = 90,  rxx = 0.80, ryy = 0.75)
#' mes_meta(r1, r2)
#'
#' \donttest{
#' # Random-effects meta-analysis with metafor
#' if (requireNamespace("metafor", quietly = TRUE)) {
#'   metafor::rma(yi = meta_df$yi, vi = meta_df$vi)
#' }
#' }
#'
#' @export
mes_meta <- function(..., study_id = NULL) {

  args <- list(...)

  # Unwrap a single list argument
  if (length(args) == 1L && is.list(args[[1L]]) &&
      !inherits(args[[1L]], "mes")) {
    args <- args[[1L]]
  }

  if (!all(vapply(args, inherits, logical(1L), what = "mes"))) {
    stop("All inputs must be of class 'mes'.", call. = FALSE)
  }

  k <- length(args)

  if (is.null(study_id)) {
    study_id <- paste0("Study_", seq_len(k))
  } else if (length(study_id) != k) {
    stop("'study_id' length must match the number of studies.",
         call. = FALSE)
  }

  rows <- vector("list", k)
  for (i in seq_len(k)) {
    obj <- args[[i]]

    yi <- obj$estimate

    # invert only symmetric CIs (HS, Bonett); BCa intervals are asymmetric,
    # so bootstrap objects take the analytical fallback below
    analytic_ci <- !is.null(obj$ci_method) && !is.na(obj$ci_method) &&
                   obj$ci_method %in% c("hunter_schmidt", "bonett")

    vi <- if (!is.null(obj$ci) && all(!is.na(obj$ci)) && analytic_ci) {
      alpha <- 1 - obj$conf_level
      crit <- if (obj$ci_method == "bonett" &&
                  !is.null(obj$n$n1) && !is.null(obj$n$n2)) {
        stats::qt(1 - alpha / 2, df = obj$n$n1 + obj$n$n2 - 2L)
      } else {
        stats::qnorm(1 - alpha / 2)
      }
      ((obj$ci[2L] - obj$ci[1L]) / (2 * crit))^2
    } else if (obj$type %in% c("d", "g") &&
               !is.null(obj$n$n1) && !is.null(obj$n$n2)) {
      # Hedges & Olkin (1985) variance; for g, Var(g*) = J^2 * Var(d) / rxx
      n1 <- obj$n$n1; n2 <- obj$n$n2
      if (obj$type == "g" && !is.null(obj$extras$cohen_d_raw) &&
          !is.null(obj$extras$hedges_J)) {
        d_underlying <- obj$extras$cohen_d_raw
        J            <- obj$extras$hedges_J
      } else {
        d_underlying <- obj$raw_estimate
        J            <- 1
      }
      var_d <- (n1 + n2) / (n1 * n2) +
               d_underlying^2 / (2 * (n1 + n2))
      var_dstar <- var_d / obj$reliability$rxx
      if (!is.null(obj$reliability$rxx_se)) {
        var_dstar <- var_dstar + d_underlying^2 *
          obj$reliability$rxx_se^2 / (4 * obj$reliability$rxx^3)
      }
      (J^2) * var_dstar
    } else if (obj$type == "r" && !is.null(obj$n$total)) {
      # Corrected correlation: large-sample Var(r) = (1 - r^2)^2 / (n - 1),
      # scaled by the fixed attenuation factor rxx * ryy. A one-sided
      # correction (no ryy) treats the second variable as error-free.
      r_obs    <- obj$raw_estimate
      n_r      <- obj$n$total
      ryy_val  <- if (!is.null(obj$reliability$ryy)) obj$reliability$ryy else 1
      atten_sq <- obj$reliability$rxx * ryy_val
      ((1 - r_obs^2)^2 / (n_r - 1)) / atten_sq
    } else NA_real_

    rxx_val <- if (!is.null(obj$reliability$rxx))
                 obj$reliability$rxx
               else if (!is.null(obj$reliability$marginal_reliability))
                 obj$reliability$marginal_reliability
               else NA_real_

    rows[[i]] <- data.frame(
      study = study_id[i],
      yi    = yi,
      vi    = vi,
      sei   = if (!is.na(vi)) sqrt(vi) else NA_real_,
      ni    = if (!is.null(obj$n$total)) obj$n$total else NA_integer_,
      rxx   = rxx_val,
      type  = obj$type
    )
  }

  do.call(rbind, rows)
}
