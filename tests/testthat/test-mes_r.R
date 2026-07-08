test_that("mes_r point estimate matches disattenuate and CI the formula", {
  fit <- mes_r(r = 0.30, n = 120, rxx = 0.85, ryy = 0.78)

  expect_s3_class(fit, "mes")
  expect_equal(fit$type, "r")
  expect_equal(fit$estimate, 0.30 / sqrt(0.85 * 0.78), tolerance = 1e-10)
  expect_equal(fit$raw_estimate, 0.30)

  z    <- atanh(0.30)
  se_z <- 1 / sqrt(120 - 3)
  crit <- qnorm(0.975)
  expected <- tanh(c(z - crit * se_z, z + crit * se_z)) / sqrt(0.85 * 0.78)
  expect_equal(fit$ci, expected, tolerance = 1e-10)
  expect_equal(fit$ci_method, "fisher_z")
})

test_that("mes_r one-sided correction warns and divides by sqrt(rxx)", {
  expect_warning(fit <- mes_r(r = 0.30, n = 100, rxx = 0.85),
                 regexp = "one-sided")
  expect_equal(fit$estimate, 0.30 / sqrt(0.85), tolerance = 1e-10)
})

test_that("mes_r truncates corrected bounds to [-1, 1]", {
  fit <- suppressWarnings(mes_r(r = 0.90, n = 30, rxx = 0.85, ryy = 0.85))
  expect_lte(fit$ci[2], 1)
  expect_gte(fit$ci[1], -1)
})

test_that("mes_r validates inputs", {
  expect_error(mes_r(r = 1.2, n = 100, rxx = 0.8), regexp = "between")
  expect_error(mes_r(r = 0.3, n = 3, rxx = 0.8), regexp = "at least 4")
  expect_error(mes_r(r = 0.3, n = 100, rxx = 1.5))
})

test_that("mes_r objects print and report without error", {
  fit <- mes_r(r = 0.40, n = 150, rxx = 0.80, ryy = 0.75)
  expect_output(print(fit))
  out <- mes_report(fit, style = "inline", language = "en")
  expect_true(grepl("corrected r", out))
})

test_that("mes_report includes both reliabilities for two-sided r", {
  fit <- mes_r(r = 0.40, n = 150, rxx = 0.80, ryy = 0.75)
  out_en <- mes_report(fit, style = "inline", language = "en")
  expect_true(grepl("r_xx = 0.80", out_en))
  expect_true(grepl("r_yy = 0.75", out_en))
  out_tr <- mes_report(fit, style = "inline", language = "tr")
  expect_true(grepl("r_yy = 0.75", out_tr))
})
