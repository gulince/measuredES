test_that("mes_ci switches CI method on existing object", {
  fit <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(80, 80),
               rxx = 0.85, ci = TRUE, ci_method = "hunter_schmidt")

  fit_bonett <- mes_ci(fit, method = "bonett")
  expect_s3_class(fit_bonett, "mes")
  expect_equal(fit_bonett$ci_method, "bonett")
  expect_length(fit_bonett$ci, 2L)
  expect_true(fit_bonett$ci[1] < fit_bonett$estimate)
  expect_true(fit_bonett$ci[2] > fit_bonett$estimate)
})

test_that("mes_ci hunter_schmidt gives reasonable interval", {
  fit <- mes_d(means = c(0, 0.5), sds = c(1, 1), ns = c(100, 100), rxx = 0.80)
  fit_hs <- mes_ci(fit, method = "hunter_schmidt")
  expect_true(all(is.finite(fit_hs$ci)))
  expect_true(fit_hs$ci[1] < fit_hs$ci[2])
})

test_that("mes_ci errors for bootstrap without raw data", {
  fit <- mes_d(means = c(100, 108), sds = c(15, 15), ns = c(60, 60), rxx = 0.85)
  expect_error(mes_ci(fit, method = "bootstrap"))
})

test_that("mes_ci inherits conf_level from the input object", {
  fit90 <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
                 rxx = 0.85, ci = TRUE, conf_level = 0.90,
                 ci_method = "hunter_schmidt")
  fit90_bon <- mes_ci(fit90, method = "bonett")
  expect_equal(fit90_bon$conf_level, 0.90)

  # explicit conf_level still overrides
  fit90_bon99 <- mes_ci(fit90, method = "bonett", conf_level = 0.99)
  expect_equal(fit90_bon99$conf_level, 0.99)
})

test_that("mes_ci on a 'g' object matches mes_g()'s own CI", {
  # mes_ci must apply J, not treat g_raw as Cohen's d
  args <- list(means = c(100, 110), sds = c(15, 15), ns = c(60, 60),
               rxx = 0.85)
  fit_g_hs     <- do.call(mes_g, c(args, list(ci = TRUE,
                                              ci_method = "hunter_schmidt")))
  fit_g_via_ci <- mes_ci(do.call(mes_g, c(args, list(ci = FALSE))),
                         method = "hunter_schmidt")
  expect_equal(fit_g_via_ci$ci, fit_g_hs$ci, tolerance = 1e-8)

  fit_g_bon     <- do.call(mes_g, c(args, list(ci = TRUE,
                                               ci_method = "bonett")))
  fit_g_via_bon <- mes_ci(do.call(mes_g, c(args, list(ci = FALSE))),
                          method = "bonett")
  expect_equal(fit_g_via_bon$ci, fit_g_bon$ci, tolerance = 1e-8)
})
