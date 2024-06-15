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
