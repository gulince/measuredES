test_that("reliability_check accepts valid values", {
  expect_invisible(reliability_check(0.85))
  expect_invisible(reliability_check(0.60))
  expect_invisible(reliability_check(0.99))
})

test_that("reliability_check errors on invalid values", {
  expect_error(reliability_check(1.5))
  expect_error(reliability_check(-0.1))
  expect_error(reliability_check(0))
})

test_that("reliability_check returns input invisibly", {
  result <- reliability_check(0.78)
  expect_equal(result, 0.78)
})
