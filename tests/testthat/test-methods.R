test_that("plot.mes runs without error", {
  fit <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
               rxx = 0.80, ci = TRUE)
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) })
  result <- plot(fit)
  # ggplot2 path returns a ggplot object; base-R path returns invisibly (NULL)
  expect_true(is.null(result) || inherits(result, "ggplot"))
})

test_that("print.mes produces output without error", {
  fit <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
               rxx = 0.80, ci = FALSE)
  expect_output(print(fit))
})

test_that("summary.mes produces output without error", {
  fit <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
               rxx = 0.80, ci = TRUE)
  expect_output(summary(fit))
})

test_that("as.data.frame.mes returns a data frame with expected columns", {
  fit <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
               rxx = 0.80, ci = TRUE)
  df <- as.data.frame(fit)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("raw_estimate", "corrected_estimate", "ci_low", "ci_high") %in%
                    names(df)))
})
