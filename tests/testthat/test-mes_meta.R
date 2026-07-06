test_that("mes_meta combines multiple mes_d objects", {
  set.seed(2026)
  studies <- lapply(1:5, function(i) {
    mes_d(means = c(0, 0.5 + rnorm(1, 0, 0.1)),
          sds   = c(1, 1),
          ns    = c(40 + 10*i, 40 + 10*i),
          rxx   = 0.80,
          ci    = TRUE,
          ci_method = "hunter_schmidt")
  })

  res <- mes_meta(studies)

  expect_s3_class(res, "data.frame")
  expect_true(all(c("yi", "vi") %in% names(res)))
  expect_equal(nrow(res), 5L)
  expect_true(all(res$vi > 0))
})

test_that("mes_meta errors on non-mes input", {
  expect_error(mes_meta(list(1, 2, 3)))
  expect_error(mes_meta("not a list"))
})

test_that("mes_meta falls back to analytical variance when the CI is NA", {
  fit <- mes_d(means = c(100, 110), sds = c(15, 15), ns = c(50, 50),
               rxx = 0.85, ci = FALSE)
  fit$ci <- c(NA_real_, NA_real_)   # simulate a failed CI
  res <- mes_meta(fit)
  expect_true(is.finite(res$vi))
  expect_true(res$vi > 0)
})

test_that("mes_meta no-CI fallback uses correct variance for 'g' objects", {
  # fallback must use Var(g*) = J^2 * Var(d) / rxx, not d's variance at g_raw
  args <- list(means = c(100, 110), sds = c(15, 15), ns = c(80, 80),
               rxx = 0.85)
  g_with_ci <- do.call(mes_g, c(args, list(ci = TRUE,
                                           ci_method = "hunter_schmidt")))
  g_no_ci   <- do.call(mes_g, c(args, list(ci = FALSE)))

  vi_with <- mes_meta(g_with_ci)$vi
  vi_no   <- mes_meta(g_no_ci)$vi

  expect_equal(vi_with, vi_no, tolerance = 1e-8)
})

test_that("mes_meta falls through to analytical variance for bootstrap CIs", {
  # BCa is asymmetric, so vi must come from the analytical fallback (= no-CI)
  set.seed(2026)
  x1 <- rnorm(40, 0, 1); x2 <- rnorm(40, 0.4, 1)
  fit_boot <- mes_d(x1, x2, rxx = 0.85,
                    ci_method = "bootstrap", n_boot = 500)
  fit_noci <- mes_d(x1, x2, rxx = 0.85, ci = FALSE)

  vi_boot <- mes_meta(fit_boot)$vi
  vi_noci <- mes_meta(fit_noci)$vi

  expect_equal(vi_boot, vi_noci, tolerance = 1e-8)
})

test_that("mes_meta computes a finite variance for correlation objects", {
  r1 <- mes_r(r = 0.30, n = 120, rxx = 0.85, ryy = 0.78)
  r2 <- mes_r(r = 0.45, n = 90,  rxx = 0.80, ryy = 0.75)

  res <- mes_meta(r1, r2, study_id = c("A", "B"))

  expect_equal(nrow(res), 2L)
  expect_true(all(is.finite(res$vi)))
  expect_true(all(res$vi > 0))
  expect_true(all(is.finite(res$sei)))
  expect_equal(res$sei, sqrt(res$vi), tolerance = 1e-12)

  # closed form: Var(r*) = (1 - r^2)^2 / ((n - 1) * rxx * ryy)
  expected_vi1 <- (1 - 0.30^2)^2 / ((120 - 1) * 0.85 * 0.78)
  expect_equal(res$vi[1L], expected_vi1, tolerance = 1e-8)
})

test_that("mes_meta handles one-sided corrected correlations (ryy absent)", {
  r1  <- suppressWarnings(mes_r(r = 0.40, n = 100, rxx = 0.90))
  res <- mes_meta(r1)

  expect_true(is.finite(res$vi))
  expected_vi <- (1 - 0.40^2)^2 / ((100 - 1) * 0.90)
  expect_equal(res$vi, expected_vi, tolerance = 1e-8)
})

test_that("mes_meta correlation output is compatible with metafor", {
  skip_if_not_installed("metafor")

  r_studies <- list(
    mes_r(r = 0.28, n = 110, rxx = 0.83, ryy = 0.79),
    mes_r(r = 0.35, n = 140, rxx = 0.88, ryy = 0.81),
    mes_r(r = 0.41, n = 95,  rxx = 0.80, ryy = 0.76)
  )
  res <- mes_meta(r_studies)
  fit <- metafor::rma(yi = yi, vi = vi, data = res, method = "REML")
  expect_s3_class(fit, "rma")
  expect_true(is.finite(fit$beta))
})

test_that("mes_meta result is compatible with metafor", {
  skip_if_not_installed("metafor")

  set.seed(2026)
  studies <- lapply(1:6, function(i) {
    mes_d(means = c(0, 0.4),
          sds   = c(1, 1),
          ns    = c(50, 50),
          rxx   = 0.85,
          ci    = TRUE,
          ci_method = "hunter_schmidt")
  })

  res <- mes_meta(studies)
  fit <- metafor::rma(yi = yi, vi = vi, data = res, method = "REML")
  expect_s3_class(fit, "rma")
  expect_true(is.finite(fit$beta))
})
