# Argument validation shared by mes_d(), mes_g(), mes_r() and mes_ci()

x1 <- c(4, 5, 6, 7, 8, 5, 6)
x2 <- c(6, 7, 8, 9, 10, 8, 9)

test_that("conf_level must be a single number strictly inside (0, 1)", {
  for (bad in list(NA_real_, 0, 1, 1.2, c(0.90, 0.95), "0.95")) {
    expect_error(mes_d(x1, x2, rxx = 0.8, conf_level = bad), "conf_level")
    expect_error(mes_g(x1, x2, rxx = 0.8, conf_level = bad), "conf_level")
    expect_error(mes_r(0.3, 50, rxx = 0.8, ryy = 0.8, conf_level = bad),
                 "conf_level")
  }
  fit <- mes_d(x1, x2, rxx = 0.8, ci = FALSE)
  expect_error(mes_ci(fit, method = "bonett", conf_level = NA_real_),
               "conf_level")
})

test_that("n_boot must be a whole number of at least 10 when bootstrapping", {
  for (bad in list(NA_real_, 0, -5, 9, 10.5, c(100, 200))) {
    expect_error(mes_d(x1, x2, rxx = 0.8, ci_method = "bootstrap",
                       n_boot = bad), "n_boot")
    expect_error(mes_g(x1, x2, rxx = 0.8, ci_method = "bootstrap",
                       n_boot = bad), "n_boot")
  }
  # n_boot is not used (and so not checked) for analytical intervals
  expect_silent(mes_d(x1, x2, rxx = 0.8, n_boot = 0))
})

test_that("summary statistics must be complete, finite and plausible", {
  ok <- function(means = c(5, 8), sds = c(1, 1), ns = c(10, 10)) {
    mes_d(means = means, sds = sds, ns = ns, rxx = 0.8)
  }
  expect_s3_class(ok(), "mes")
  expect_error(ok(means = c(NA, 8)), "finite and non-missing")
  expect_error(ok(sds = c(1, NA)), "finite and non-missing")
  expect_error(ok(ns = c(10, NA)), "finite and non-missing")
  expect_error(ok(means = c(Inf, 8)), "finite and non-missing")
  expect_error(ok(sds = c(1, 0)), "positive")
  expect_error(ok(ns = c(10.5, 10)), "whole numbers")
  expect_error(ok(ns = c(1, 10)), "at least 2")
  expect_error(ok(means = 5), "length-2")
  expect_error(ok(means = c("a", "b")), "numeric")
})
