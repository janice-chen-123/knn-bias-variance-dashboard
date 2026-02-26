# Unit tests for the STA380 k-NN Monte Carlo functions
# These tests are intended to be run with testthat.
#
# In a package, these would typically live under:
#   tests/testthat/test-knn-mc.R
#
# For the Basic Implementation Check, a single test file is acceptable.

library(testthat)

# If developing as a package, you would use:
# devtools::load_all()
# Here, source the functions directly:
source("function_2.0.R")

test_that("dgp_x returns correct dimensions", {
  x1 <- dgp_x(n = 10, d = 1, seed = 1)
  x2 <- dgp_x(n = 10, d = 2, seed = 1)
  expect_equal(dim(x1), c(10, 1))
  expect_equal(dim(x2), c(10, 2))
})

test_that("dgp_f matches known special values (1d)", {
  expect_equal(dgp_f(0, type = "1d"), 0, tolerance = 1e-12)
  expect_equal(dgp_f(0.5, type = "1d"), 0, tolerance = 1e-12) # sin(pi)=0
})

test_that("knn_predict with k=n equals the training mean for any test point", {
  x_train <- matrix(c(0, 0.5, 1), ncol = 1)
  y_train <- c(1, 2, 3)
  x_test <- c(0.2, 0.8)
  pred <- knn_predict(x_train, y_train, x_test, k = length(y_train))
  expect_equal(pred, rep(mean(y_train), length(x_test)))
})

test_that("mc_knn_decomp returns expected structure and nonnegative components", {
  res <- mc_knn_decomp(n = 30, k = 3, B = 40, sigma = 0.2, d = 1, type = "1d", grid_size = 21, seed = 123)
  expect_true(is.list(res))
  expect_true(all(c("settings", "grid_axis", "pointwise") %in% names(res)))

  pw <- res$pointwise
  expect_equal(nrow(pw), 21)
  expect_true(all(pw$variance >= 0))
  expect_true(all(pw$bias2 >= 0))
  expect_true(all(abs(pw$mse_f - (pw$bias2 + pw$variance)) < 1e-10))
  expect_true(all(abs(pw$mse_y - (pw$mse_f + res$settings$sigma^2)) < 1e-10))
})

test_that("mc_knn_curve returns one row per k", {
  curve <- mc_knn_curve(k_values = c(1, 3, 5), n = 30, B = 20, sigma = 0.2, d = 1, type = "1d", grid_size = 21, seed = 1)
  expect_equal(nrow(curve), 3)
  expect_true(all(curve$k %in% c(1, 3, 5)))
})
