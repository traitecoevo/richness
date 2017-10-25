context("Zipf")

n <- floor(runif(20,1,10)^2)
f1 <-  fitZipf(n)
vars <- c("distribution","exponent","fit")

test_that("Fit Zipf distribution", {
  expect_named(f1)
  expect_is(f1, "list")
  expect_length(f1, 3L)
  expect_identical(names(f1), vars)
  expect_length(f1$distribution, length(n))

  expect_error(fitZipf(c(n, NA)))
  expect_error(fitZipf(NULL))

})
