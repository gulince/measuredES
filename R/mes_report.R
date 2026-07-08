# mes_report.R: APA-7 text and table output from a mes object


#' APA-7 style reporting for 'mes' objects
#'
#' @description
#' Converts an 'mes' object to text following the recommendations of the
#' APA Publication Manual (7th edition). Produces both in-line text and
#' tabular formats.
#'
#' @details
#' APA style calls for effect sizes to be reported and, where possible,
#' accompanied by confidence intervals (American Psychological Association,
#' 2020). In line with long-standing transparency guidelines (Wilkinson &
#' APA Task Force on Statistical Inference, 1999), the inline output always
#' shows the observed estimate next to the corrected one, so readers can
#' see that (and how much) a correction was applied.
#'
#' @param x An 'mes' object.
#' @param style Output style: \code{"inline"} (a sentence-style string),
#'   \code{"table"} (a data frame), or \code{"both"} (a list with both).
#' @param language Output language: \code{"en"} (English) or \code{"tr"}
#'   (Turkish).
#' @param digits Number of decimal places.
#'
#' @return Character text, data frame, or list.
#'
#' @references
#' American Psychological Association. (2020). \emph{Publication manual of
#' the American Psychological Association} (7th ed.).
#'
#' Wilkinson, L., & APA Task Force on Statistical Inference. (1999).
#' Statistical methods in psychology journals: Guidelines and explanations.
#' \emph{American Psychologist}, 54(8), 594--604.
#'
#' @examples
#' set.seed(2026)
#' a <- rnorm(50, 100, 15); b <- rnorm(50, 110, 15)
#' fit <- mes_d(a, b, rxx = 0.82)
#'
#' mes_report(fit, style = "inline", language = "en")
#' mes_report(fit, style = "inline", language = "tr")
#' mes_report(fit, style = "table")
#'
#' @export
mes_report <- function(x,
                       style = c("inline", "table", "both"),
                       language = c("en", "tr"),
                       digits = 2L) {

  if (!inherits(x, "mes")) {
    stop("'x' must be an object of class 'mes'.", call. = FALSE)
  }

  style    <- match.arg(style)
  language <- match.arg(language)

  symbol <- switch(x$type,
                   d = "d",
                   g = "g",
                   r = "r",
                   "ES")

  fmt <- function(v) formatC(v, digits = digits, format = "f")

  inline_text <- NULL
  if (style %in% c("inline", "both")) {
    if (language == "en") {
      inline_text <- sprintf(
        "corrected %s = %s (observed %s = %s)",
        symbol, fmt(x$estimate), symbol, fmt(x$raw_estimate))

      if (!is.null(x$ci) && all(!is.na(x$ci))) {
        inline_text <- paste0(inline_text,
                              sprintf(", %d%% CI [%s, %s]",
                                      round(x$conf_level * 100),
                                      fmt(x$ci[1L]), fmt(x$ci[2L])))
      }

      if (!is.null(x$reliability$rxx)) {
        inline_text <- paste0(inline_text,
                              sprintf("; reliability r_xx = %s",
                                      fmt(x$reliability$rxx)))
        if (!is.null(x$reliability$ryy)) {
          inline_text <- paste0(inline_text,
                                sprintf(", r_yy = %s",
                                        fmt(x$reliability$ryy)))
        }
      } else if (!is.null(x$reliability$marginal_reliability)) {
        inline_text <- paste0(inline_text,
                              sprintf("; marginal reliability = %s",
                                      fmt(x$reliability$marginal_reliability)))
      }
    } else {
      inline_text <- sprintf(
        "D\u00fczeltilmi\u015f %s = %s (g\u00f6zlenen %s = %s)",
        symbol, fmt(x$estimate), symbol, fmt(x$raw_estimate))

      if (!is.null(x$ci) && all(!is.na(x$ci))) {
        inline_text <- paste0(inline_text,
                              sprintf(", %%%d GA [%s, %s]",
                                      round(x$conf_level * 100),
                                      fmt(x$ci[1L]), fmt(x$ci[2L])))
      }
      if (!is.null(x$reliability$rxx)) {
        inline_text <- paste0(inline_text,
                              sprintf("; g\u00fcvenirlik r_xx = %s",
                                      fmt(x$reliability$rxx)))
        if (!is.null(x$reliability$ryy)) {
          inline_text <- paste0(inline_text,
                                sprintf(", r_yy = %s",
                                        fmt(x$reliability$ryy)))
        }
      }
    }
  }

  table_out <- NULL
  if (style %in% c("table", "both")) {
    table_out <- as.data.frame(x)
  }

  switch(style,
         inline = inline_text,
         table  = table_out,
         both   = list(inline = inline_text, table = table_out))
}
