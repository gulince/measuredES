test_that("mes_g works with raw data and applies Hedges correction", {
  set.seed(2026)
  x1 <- rnorm(40, 100, 15)
  x2 <- rnorm(40, 110, 15)
  fit <- mes_g(x1 = x1, x2 = x2, rxx = 0.85, ci = FALSE)

  expect_s3_class(fit, "mes")
  expect_equal(fit$type, "g")
  expect_true(fit$estimate < 0)                                # x1 < x2
  expect_true(abs(fit$estimate) > abs(fit$raw_estimate))        # correction inflates
  # |g| < |d| because J < 1
  fit_d <- mes_d(x1 = x1, x2 = x2, rxx = 0.85, ci = FALSE)
  expect_true(abs(fit$raw_estimate) < abs(fit_d$raw_estimate))
})

test_that("mes_g works with summary statistics", {
  fit <- mes_g(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
               rxx = 0.80, ci = FALSE)
  expect_s3_class(fit, "mes")
  # raw_estimate = J * (-10/15)
  df <- 50 + 50 - 2
  J  <- 1 - 3 / (4 * df - 1)
  expect_equal(fit$raw_estimate, J * (-10/15), tolerance = 1e-6)
  expect_equal(fit$estimate, fit$raw_estimate / sqrt(0.80), tolerance = 1e-6)
})

test_that("mes_g errors on invalid reliability values", {
  expect_error(mes_g(means = c(0, 1), sds = c(1, 1), ns = c(10, 10),
                     rxx = 1.5))
  expect_error(mes_g(means = c(0, 1), sds = c(1, 1), ns = c(10, 10),
                     rxx = -0.1))
})

test_that("mes_g returns an object with valid CI when requested", {
  set.seed(2026)
  fit <- mes_g(means = c(100, 110), sds = c(15, 15), ns = c(80, 80),
               rxx = 0.85, ci = TRUE, ci_method = "hunter_schmidt")
  expect_length(fit$ci, 2L)
  expect_true(fit$ci[1] < fit$estimate)
  expect_true(fit$ci[2] > fit$estimate)
})

test_that("mes_g Bonett CI produces valid bounds", {
  fit <- mes_g(means = c(100, 110), sds = c(15, 15), ns = c(80, 80),
               rxx = 0.85, ci = TRUE, ci_method = "bonett")
  expect_length(fit$ci, 2L)
  expect_true(fit$ci[1] < fit$estimate)
  expect_true(fit$ci[2] > fit$estimate)
  expect_equal(fit$ci_method, "bonett")
})

test_that("mes_g Bootstrap CI produces valid bounds with raw data", {
  set.seed(2026)
  x1 <- rnorm(60, 100, 15)
  x2 <- rnorm(60, 110, 15)
  fit <- mes_g(x1 = x1, x2 = x2, rxx = 0.85, ci = TRUE,
               ci_method = "bootstrap", n_boot = 1000)
  expect_length(fit$ci, 2L)
  expect_true(fit$ci[1] < fit$estimate)
  expect_true(fit$ci[2] > fit$estimate)
  expect_equal(fit$ci_method, "bootstrap")
})

test_that("mes_g Bootstrap CI warns and falls back without raw data", {
  expect_warning(
    mes_g(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
          rxx = 0.85, ci = TRUE, ci_method = "bootstrap"),
    regexp = "raw data"
  )
})

test_that("mes_g bootstrap drops NAs before resampling", {
  # same seed on pre-cleaned data must give the identical CI
  base1 <- rnorm(40, 100, 15)
  base2 <- rnorm(40, 110, 15)

  set.seed(99)
  fit_na    <- mes_g(c(base1, NA, NA), c(base2, NA), rxx = 0.85,
                     ci_method = "bootstrap", n_boot = 400)
  set.seed(99)
  fit_clean <- mes_g(base1, base2, rxx = 0.85,
                     ci_method = "bootstrap", n_boot = 400)

  expect_equal(fit_na$n$n1, 40L)
  expect_equal(fit_na$n$n2, 40L)
  expect_equal(fit_na$ci, fit_clean$ci, tolerance = 1e-12)
})

test_that("mes_g supports the formula interface", {
  set.seed(2026)
  df <- data.frame(
    score = c(rnorm(40, 100, 15), rnorm(40, 110, 15)),
    group = rep(c("control", "treatment"), each = 40)
  )
  fit_formula <- mes_g(formula = score ~ group, data = df, rxx = 0.85)
  fit_raw     <- mes_g(x1 = df$score[df$group == "control"],
                       x2 = df$score[df$group == "treatment"],
                       rxx = 0.85)
  expect_s3_class(fit_formula, "mes")
  expect_equal(fit_formula$estimate, fit_raw$estimate, tolerance = 1e-10)
})
