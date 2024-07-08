test_that("check ddnorm", {
    #sums to 1
    expect_equal(sum(ddnorm(-10:10)), 1)

    #symmetric distribution
    expect_equal(ddnorm(-c(1:10)), ddnorm(1:10))

    #center moves with mean
    expect_equal(ddnorm(1, 0, 1), ddnorm(4, 3, 1))
    
    expect_true(memoise::is.memoised(ddnorm_constant))

})

test_that("check rddnorm", {
    set.seed(1)
    n <- 10000
    mu <- 300
    sigma <- 2
    smpl <- rdnorm(n, mu, sigma)

    #check mean is mu
    expect_equal(mean(smpl), 300, tolerance = 1e-2)

    #check sample are integer valued
    expect_equal(sum(abs(ceiling(smpl) - floor(smpl))), 0)
})
