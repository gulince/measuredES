# mes_ci.R: recompute the CI on an existing mes object


#' Confidence interval for a corrected effect size
#'
#' @description
#' Swap the CI method on an existing \code{mes} object without recomputing
#' the effect size itself. Useful when a reviewer asks for a different
#' interval (e.g., switching from Hunter-Schmidt to Bonett).
#'
#' @param object An 'mes' object.
#' @param method CI method: \code{"hunter_schmidt"}, \code{"bonett"}, or
#'   \code{"bootstrap"}.
#' @param conf_level Confidence level. Defaults to the level stored on the
#'   incoming object, so the method change alone doesn't quietly change the
#'   coverage.
#' @param n_boot Number of bootstrap replications (only for bootstrap).
#'
#' @return Updated 'mes' object.
#'
#' @examples
#' set.seed(2026)
#' a <- rnorm(40, 100, 15); b <- rnorm(40, 108, 15)
#' fit <- mes_d(a, b, rxx = 0.85, ci = FALSE)
#' fit_ci <- mes_ci(fit, method = "bonett")
#' fit_ci$ci
#'
#' @export
mes_ci <- function(object,
                   method = c("hunter_schmidt", "bonett", "bootstrap"),
                   conf_level = NULL,
                   n_boot = 2000L) {

  if (!inherits(object, "mes")) {
    stop("'object' must be of class 'mes'.", call. = FALSE)
  }

  method <- match.arg(method)

  if (is.null(conf_level)) {
    conf_level <- object$conf_level
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be in (0, 1).", call. = FALSE)
  }

  if (!object$type %in% c("d", "g")) {
    stop("mes_ci only works with 'd' and 'g' objects for now. ",
         "For other types, pass ci = TRUE directly to the original function.",
         call. = FALSE)
  }

  n1     <- object$n$n1
  n2     <- object$n$n2
  rxx    <- object$reliability$rxx
  rxx_se <- object$reliability$rxx_se
  s1     <- object$extras$s1
  s2     <- object$extras$s2

  # for g objects: Var(g*) = J^2 * Var(d) / rxx, with d_raw and J from extras
  if (object$type == "g") {
    d_raw <- object$extras$cohen_d_raw
    J     <- object$extras$hedges_J
    if (is.null(d_raw) || is.null(J)) {
      stop("This 'g' object is missing the d_raw/J information needed to ",
           "recompute its CI. Re-run mes_g() instead.", call. = FALSE)
    }
  } else {
    d_raw <- object$raw_estimate
    J     <- 1
  }

  d_ci <- if (method == "hunter_schmidt") {
    .hs_ci_d(d_raw, rxx, n1, n2, conf_level, rxx_se)
  } else if (method == "bonett") {
    if (is.null(s1) || is.null(s2)) {
      stop("Bonett's variance needs the two group SDs (s1, s2). The 'mes' ",
           "object doesn't carry them; recompute with mes_d/mes_g.",
           call. = FALSE)
    }
    .bonett_ci_d(d_raw, rxx, n1, n2, s1, s2, conf_level, rxx_se)
  } else {
    stop("Bootstrap CI requires raw data; use mes_d(..., ci_method = ",
         "'bootstrap') or mes_g(..., ci_method = 'bootstrap') instead.",
         call. = FALSE)
  }

  # g = J * d, so the g interval is J times the d interval
  object$ci         <- J * d_ci
  object$ci_method  <- method
  object$conf_level <- conf_level

  object
}
