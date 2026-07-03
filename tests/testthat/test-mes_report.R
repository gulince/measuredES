test_that("mes_report produces inline English text", {
  fit <- mes_d(means = c(100, 108), sds = c(15, 15), ns = c(60, 60),
               rxx = 0.85, ci = TRUE, ci_method = "hunter_schmidt")
  out <- mes_report(fit, style = "inline", language = "en")
  expect_type(out, "character")
  expect_true(grepl("corrected d", out))
  expect_true(grepl("0.85", out))
})

test_that("mes_report produces inline Turkish text", {
  fit <- mes_d(means = c(100, 108), sds = c(15, 15), ns = c(60, 60),
               rxx = 0.80, ci = FALSE)
  out <- mes_report(fit, style = "inline", language = "tr")
  expect_type(out, "character")
  expect_true(grepl("zeltilmi", out, ignore.case = TRUE) ||
              grepl("d =", out))
})

test_that("mes_report table style returns data frame", {
  fit <- mes_d(means = c(0, 0.5), sds = c(1, 1), ns = c(50, 50), rxx = 0.80)
  out <- mes_report(fit, style = "table")
  expect_true(is.data.frame(out) || is.list(out))
})
