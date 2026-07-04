# methods.R: print, summary, plot and as.data.frame for the mes class


#' Print method for 'mes' objects
#'
#' @param x An object of class 'mes'.
#' @param digits Number of decimal places for display.
#' @param ... Further arguments (unused).
#'
#' @export
print.mes <- function(x, digits = 3, ...) {

  type_label <- switch(x$type,
                       "d"      = "Cohen's d (corrected)",
                       "g"      = "Hedges' g (corrected)",
                       "r"      = "Pearson r (corrected)",
                       "Effect size")

  cat("\n --", type_label, "--\n\n")

  cat(sprintf(" Raw (observed) estimate  : %.*f\n",
              digits, x$raw_estimate))
  cat(sprintf(" Corrected estimate       : %.*f\n",
              digits, x$estimate))

  if (!is.null(x$ci) && all(!is.na(x$ci))) {
    cat(sprintf(" %d%% Confidence interval  : [%.*f, %.*f]\n",
                round(x$conf_level * 100),
                digits, x$ci[1L],
                digits, x$ci[2L]))
    cat(sprintf(" CI method                : %s\n", x$ci_method))
  }

  if (length(x$reliability) > 0L) {
    cat("\n Reliability information:\n")
    for (nm in names(x$reliability)) {
      cat(sprintf("   %s = %.*f\n", nm, digits, x$reliability[[nm]]))
    }
  }

  if (length(x$n) > 0L) {
    cat("\n Sample size information:\n")
    for (nm in names(x$n)) {
      cat(sprintf("   %s = %d\n", nm, x$n[[nm]]))
    }
  }

  if (!is.na(x$interpretation)) {
    cat(sprintf("\n Cohen (1988) benchmark: %s\n", x$interpretation))
  }

  cat("\n")
  invisible(x)
}


#' Summary method for 'mes' objects
#'
#' @param object An object of class 'mes'.
#' @param ... Further arguments.
#'
#' @export
summary.mes <- function(object, ...) {
  print(object, ...)

  if (!is.na(object$raw_estimate) && object$raw_estimate != 0) {
    inflation <- (object$estimate - object$raw_estimate) / object$raw_estimate
    cat(sprintf(" Correction changed the estimate by %.1f%%\n\n",
                inflation * 100))
  }

  invisible(object)
}


#' Convert an 'mes' object to a data frame
#'
#' @description
#' Turns a \code{mes} object into a one-row data frame: raw and corrected
#' estimate, CI bounds, method, reliability, sample size, and Cohen label.
#' For a \code{metafor}-ready frame (yi, vi, sei, ni), use \code{mes_meta()}.
#'
#' @param x An 'mes' object.
#' @param row.names Row names (ignored).
#' @param optional Logical (ignored).
#' @param ... Further arguments.
#'
#' @export
as.data.frame.mes <- function(x, row.names = NULL, optional = FALSE, ...) {

  ci_low  <- if (!is.null(x$ci)) x$ci[1L] else NA_real_
  ci_high <- if (!is.null(x$ci)) x$ci[2L] else NA_real_
  n_total <- if ("total" %in% names(x$n)) x$n$total
             else if (length(x$n) > 0L) sum(unlist(x$n))
             else NA_integer_

  data.frame(
    type               = x$type,
    raw_estimate       = x$raw_estimate,
    corrected_estimate = x$estimate,
    ci_low             = ci_low,
    ci_high            = ci_high,
    ci_method          = x$ci_method,
    conf_level         = x$conf_level,
    reliability_dv     = if ("rxx" %in% names(x$reliability))
                           x$reliability$rxx else NA_real_,
    n_total            = n_total,
    interpretation     = x$interpretation
  )
}


#' Plot method for 'mes' objects
#'
#' @description
#' Displays the corrected and uncorrected effect sizes together with the
#' confidence interval in a forest-style plot. Uses \code{ggplot2} if
#' installed; otherwise falls back to base graphics.
#'
#' @param x An 'mes' object.
#' @param ... Further arguments.
#'
#' @export
plot.mes <- function(x, ...) {

  use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  if (use_ggplot) {
    df <- data.frame(
      label    = c("Raw (observed)", "Corrected"),
      estimate = c(x$raw_estimate, x$estimate),
      lower    = c(NA_real_, if (!is.null(x$ci)) x$ci[1L] else NA_real_),
      upper    = c(NA_real_, if (!is.null(x$ci)) x$ci[2L] else NA_real_)
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate,
                                          y = .data$label)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower,
                                          xmax = .data$upper),
                             width = 0.15, na.rm = TRUE,
                             orientation = "y") +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(x = "Effect size", y = NULL,
                    title = "Corrected vs. observed effect size") +
      ggplot2::theme_minimal()

    return(p)
  }

  # Base R fallback
  vals <- c(x$raw_estimate, x$estimate)
  graphics::dotchart(vals, labels = c("Raw", "Corrected"),
                     xlab = "Effect size",
                     main = "Corrected vs. observed effect size")
  graphics::abline(v = 0, lty = 2, col = "grey50")

  invisible(NULL)
}
