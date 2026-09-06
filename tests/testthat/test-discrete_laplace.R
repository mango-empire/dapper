test_that("check ddlaplace", {
  #sums to 1
  expect_equal(sum(ddlaplace(-20:20)), 1)

  #symmetric distribution
  expect_equal(ddlaplace(-c(1:20)), ddlaplace(1:20))
})

test_that("check rdlaplace", {
  set.seed(1)
  n <- 20000
  scale <- 3
  smpl <- rdlaplace(n, scale)

  #check mean is 0
  expect_equal(mean(smpl), 0, tolerance = 1e-2)

  #check sample are integer valued
  expect_equal(sum(abs(ceiling(smpl) - floor(smpl))), 0)
})

test_that("rdlaplace requires a positive integer scale", {
  for (scale in c(-1, 0, 1.5, 1 + 1e-10, Inf, NA_real_)) {
    expect_error(rdlaplace(1, scale = scale), "scale")
  }

  set.seed(1)
  for (scale in c(1, 2)) {
    draws <- rdlaplace(10, scale = scale)
    expect_length(draws, 10)
    expect_true(all(is.finite(draws) & draws == floor(draws)))
  }
})

test_that("ddlaplace supports positive noninteger scales", {
  expect_equal(sum(ddlaplace(-100:100, scale = 1.5)), 1)
})
