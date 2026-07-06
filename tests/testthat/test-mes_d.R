test_that("mes_d works correctly with raw data", {
  set.seed(2026)
  x1 <- rnorm(50, 100, 15)
  x2 <- rnorm(50, 110, 15)
  fit <- mes_d(x1 = x1, x2 = x2, rxx = 0.85, ci = FALSE)

  expect_s3_class(fit, "mes")
  expect_equal(fit$type, "d")
  expect_true(fit$estimate < 0)                              # x1 < x2
  expect_true(abs(fit$estimate) > abs(fit$raw_estimate))      # correction inflates
  expect_equal(fit$reliability$rxx, 0.85)
})

test_that("mes_d works correctly with summary statistics", {
  fit <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
               rxx = 0.80, ci = FALSE)
  expect_s3_class(fit, "mes")
  # d_observed = -10 / 15 = -0.6667
  expect_equal(fit$raw_estimate, -10/15, tolerance = 1e-6)
  # d_corrected = -0.6667 / sqrt(0.80)
  expect_equal(fit$estimate, (-10/15) / sqrt(0.80), tolerance = 1e-6)
})

test_that("mes_d errors on invalid reliability values", {
  expect_error(mes_d(means = c(0, 1), sds = c(1, 1), ns = c(10, 10),
                     rxx = 1.5))
  expect_error(mes_d(means = c(0, 1), sds = c(1, 1), ns = c(10, 10),
                     rxx = -0.1))
  expect_error(mes_d(means = c(0, 1), sds = c(1, 1), ns = c(10, 10),
                     rxx = 0))
})

test_that("Hunter-Schmidt CI produces reasonable bounds for d", {
  set.seed(2026)
  fit <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(100, 100),
               rxx = 0.80,
               ci = TRUE, ci_method = "hunter_schmidt")
  expect_length(fit$ci, 2L)
  expect_true(fit$ci[1L] < fit$estimate)
  expect_true(fit$ci[2L] > fit$estimate)
})

test_that("disattenuate general engine works correctly", {
  expect_equal(disattenuate(0.5, rxx = 0.81, type = "d"),
               0.5 / sqrt(0.81), tolerance = 1e-6)
  expect_equal(disattenuate(0.30, rxx = 0.85, ryy = 0.78, type = "r"),
               0.30 / sqrt(0.85 * 0.78), tolerance = 1e-6)
})

test_that("disattenuate(r) without ryy warns about the one-sided assumption", {
  expect_warning(disattenuate(0.30, rxx = 0.85, type = "r"),
                 regexp = "one-sided")
})

test_that("disattenuate no longer accepts eta2 / omega2", {
  expect_error(disattenuate(0.10, rxx = 0.80, type = "eta2"))
  expect_error(disattenuate(0.10, rxx = 0.80, type = "omega2"))
})

test_that("disattenuate is vectorized over 'effect' for type = 'd'", {
  d_vec <- c(0.20, 0.45, 0.63, 0.80)
  out   <- disattenuate(d_vec, rxx = 0.75, type = "d")
  expect_length(out, length(d_vec))
  expect_equal(out, d_vec / sqrt(0.75), tolerance = 1e-10)
})

test_that("disattenuate passes NA elements through unchanged", {
  d_vec <- c(0.20, NA, 0.63)
  out   <- disattenuate(d_vec, rxx = 0.80, type = "d")
  expect_equal(out, d_vec / sqrt(0.80))
  expect_true(is.na(out[2L]))
})

test_that("disattenuate is vectorized for type = 'r'", {
  r_vec <- c(0.30, 0.50)
  out   <- disattenuate(r_vec, rxx = 0.85, ryy = 0.78, type = "r")
  expect_length(out, length(r_vec))
  expect_equal(out, r_vec / sqrt(0.85 * 0.78), tolerance = 1e-10)
})

test_that("disattenuate warns once, with a count, when any |r*| > 1", {
  # 0.95 / (0.85) = 1.12 is out of bounds; 0.20 / 0.85 = 0.24 is fine
  expect_warning(
    disattenuate(c(0.20, 0.95), rxx = 0.85, ryy = 0.85, type = "r"),
    regexp = "1 corrected correlation"
  )
})

test_that("mes_d formula interface works", {
  set.seed(2026)
  df <- data.frame(
    score = c(rnorm(40, 100, 15), rnorm(40, 110, 15)),
    group = rep(c("control", "treatment"), each = 40)
  )
  fit_formula <- mes_d(formula = score ~ group, data = df, rxx = 0.85)
  fit_raw     <- mes_d(x1 = df$score[df$group == "control"],
                       x2 = df$score[df$group == "treatment"],
                       rxx = 0.85)
  expect_s3_class(fit_formula, "mes")
  expect_equal(fit_formula$estimate, fit_raw$estimate, tolerance = 1e-10)
})

test_that("mes_d formula errors on non-binary group", {
  df <- data.frame(score = rnorm(60), group = rep(letters[1:3], 20))
  expect_error(mes_d(formula = score ~ group, data = df, rxx = 0.80))
})

test_that("mes_d handles zero effect size correctly", {
  fit <- mes_d(means = c(100, 100), sds = c(15, 15), ns = c(50, 50),
               rxx = 0.80, ci = FALSE)
  expect_equal(fit$raw_estimate, 0, tolerance = 1e-10)
  expect_equal(fit$estimate, 0, tolerance = 1e-10)
})

test_that("mes_d high reliability (0.99) inflates estimate minimally", {
  fit_hi  <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
                   rxx = 0.99, ci = FALSE)
  fit_lo  <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
                   rxx = 0.50, ci = FALSE)
  expect_true(abs(fit_hi$estimate - fit_hi$raw_estimate) <
                abs(fit_lo$estimate - fit_lo$raw_estimate))
})

test_that("mes_d very low reliability warns but does not error", {
  expect_warning(
    mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
          rxx = 0.05, ci = FALSE),
    regexp = "low"
  )
})

test_that("mes_d large sample size produces narrow CI", {
  fit_small <- mes_d(means = c(100, 105), sds = c(10, 10), ns = c(30, 30),
                     rxx = 0.80, ci = TRUE, ci_method = "hunter_schmidt")
  fit_large <- mes_d(means = c(100, 105), sds = c(10, 10), ns = c(500, 500),
                     rxx = 0.80, ci = TRUE, ci_method = "hunter_schmidt")
  width_small <- fit_small$ci[2] - fit_small$ci[1]
  width_large <- fit_large$ci[2] - fit_large$ci[1]
  expect_true(width_large < width_small)
})

test_that("mes_d drops NA values from raw vectors", {
  set.seed(2026)
  x1 <- c(rnorm(30, 100, 15), NA, NA)
  x2 <- c(rnorm(30, 110, 15), NA)
  fit <- mes_d(x1 = x1, x2 = x2, rxx = 0.85, ci = FALSE)
  # n1, n2 should reflect the non-missing counts only
  expect_equal(fit$n$n1, 30L)
  expect_equal(fit$n$n2, 30L)
  expect_true(is.finite(fit$estimate))
})

test_that("mes_d errors when, after NA removal, a group has fewer than 2 obs", {
  expect_error(
    mes_d(x1 = c(1, NA, NA), x2 = c(5, 6, 7), rxx = 0.85),
    regexp = "at least"
  )
})

test_that("mes_d errors on non-numeric summary stats", {
  expect_error(
    mes_d(means = c("a", "b"), sds = c(1, 1), ns = c(20, 20), rxx = 0.85)
  )
  expect_error(
    mes_d(means = c(0, 1), sds = c(1, 1), ns = c("50", "50"), rxx = 0.85),
    regexp = "numeric"
  )
})

test_that("mes_d formula rejects more than one RHS variable", {
  df <- data.frame(y = rnorm(40),
                   g = rep(c("a", "b"), 20),
                   z = rep(c("c", "d"), each = 20))
  expect_error(
    mes_d(formula = y ~ g + z, data = df, rxx = 0.80),
    regexp = "exactly one grouping variable"
  )
})

test_that("disattenuate warns when ryy is supplied with type = 'd'", {
  expect_warning(
    disattenuate(0.5, rxx = 0.80, ryy = 0.70, type = "d"),
    regexp = "ignored"
  )
})
