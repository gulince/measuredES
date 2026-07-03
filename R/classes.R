# classes.R: the mes S3 class and its constructor


#' Internal constructor for 'mes' objects
#'
#' @description
#' Every \code{mes_*()} function packages its result through this constructor
#' so the print, summary, plot, and as.data.frame methods all behave the
#' same way.
#'
#' @param estimate Corrected effect-size estimate.
#' @param raw_estimate Uncorrected (observed) effect size.
#' @param ci Confidence interval (length-2 vector or NULL).
#' @param ci_method CI computation method.
#' @param conf_level Confidence level.
#' @param type Effect-size type: \code{"d"}, \code{"g"}, or \code{"r"}.
#' @param reliability Reliability information used (named list).
#' @param n Sample-size information (named list).
#' @param call The matched function call.
#' @param extras Optional additional fields (named list).
#'
#' @return An object of class 'mes' (a list).
#'
#' @keywords internal
.new_mes <- function(estimate,
                     raw_estimate = NA_real_,
                     ci = NULL,
                     ci_method = NA_character_,
                     conf_level = 0.95,
                     type = NA_character_,
                     reliability = list(),
                     n = list(),
                     call = NULL,
                     extras = list()) {

  obj <- list(
    estimate     = estimate,
    raw_estimate = raw_estimate,
    ci           = ci,
    ci_method    = ci_method,
    conf_level   = conf_level,
    type         = type,
    reliability  = reliability,
    n            = n,
    interpretation = .cohen_interpretation(estimate, type = type),
    call         = call,
    extras       = extras
  )

  class(obj) <- c("mes", "list")
  obj
}
