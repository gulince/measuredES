test_that("rxx_se widens the analytical interval and matches the formula", {
  fit0 <- mes_d(means = c(1, 0), sds = c(1, 1), ns = c(50, 50),
                rxx = 0.80, ci_method = "hunter_schmidt")
  fit1 <- mes_d(means = c(1, 0), sds = c(1, 1), ns = c(50, 50),
                rxx = 0.80, rxx_se = 0.05, ci_method = "hunter_schmidt")

  expect_gt(diff(fit1$ci), diff(fit0$ci))

  # hand computation of the delta-method variance
  d_raw <- fit1$raw_estimate
  var_d <- (50 + 50) / (50 * 50) + d_raw^2 / (2 * 100)
  var_c <- var_d / 0.80 + d_raw^2 * 0.05^2 / (4 * 0.80^3)
  z     <- qnorm(0.975)
  expect_equal(as.numeric(diff(fit1$ci)), 2 * z * sqrt(var_c),
               tolerance = 1e-10)
  expect_equal(fit1$reliability$rxx_se, 0.05)
})

test_that("rxx_se works with the Bonett interval too", {
  fit0 <- mes_d(means = c(1, 0), sds = c(1, 1.2), ns = c(40, 60),
                rxx = 0.75, ci_method = "bonett")
  fit1 <- mes_d(means = c(1, 0), sds = c(1, 1.2), ns = c(40, 60),
                rxx = 0.75, rxx_se = 0.06, ci_method = "bonett")
  expect_gt(diff(fit1$ci), diff(fit0$ci))
})

test_that("mes_g propagates rxx_se and mes_ci reproduces it", {
  args <- list(means = c(100, 110), sds = c(15, 15), ns = c(60, 60),
               rxx = 0.85, rxx_se = 0.04)
  g_direct <- do.call(mes_g, c(args, list(ci_method = "hunter_schmidt")))
  g_via_ci <- mes_ci(do.call(mes_g, c(args, list(ci = FALSE))),
                     method = "hunter_schmidt")
  expect_equal(g_via_ci$ci, g_direct$ci, tolerance = 1e-10)

  g_plain <- mes_g(means = c(100, 110), sds = c(15, 15), ns = c(60, 60),
                   rxx = 0.85, ci_method = "hunter_schmidt")
  expect_gt(diff(g_direct$ci), diff(g_plain$ci))
})

test_that("mes_meta stays internally consistent under rxx_se", {
  fit_ci   <- mes_d(means = c(0, 0.5), sds = c(1, 1), ns = c(80, 80),
                    rxx = 0.80, rxx_se = 0.05,
                    ci_method = "hunter_schmidt")
  fit_noci <- mes_d(means = c(0, 0.5), sds = c(1, 1), ns = c(80, 80),
                    rxx = 0.80, rxx_se = 0.05, ci = FALSE)
  expect_equal(mes_meta(fit_ci)$vi, mes_meta(fit_noci)$vi,
               tolerance = 1e-10)
})

test_that("rxx_se is rejected when invalid and ignored for bootstrap", {
  expect_error(
    mes_d(means = c(0, 1), sds = c(1, 1), ns = c(30, 30),
          rxx = 0.8, rxx_se = -0.1),
    regexp = "positive"
  )
  set.seed(2026)
  x1 <- rnorm(40); x2 <- rnorm(40, 0.4)
  expect_warning(
    mes_d(x1, x2, rxx = 0.8, rxx_se = 0.05,
          ci_method = "bootstrap", n_boot = 200),
    regexp = "ignored"
  )
})
