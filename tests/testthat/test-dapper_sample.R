test_that("check inputs", {
    dmod <- list()
    expect_error(dapper_sample(data_model = dmod))

    #sdp must be matrix or vector
    dmod <- structure(list(), class = "privacy")
    expect_error(dapper_sample(data_model = dmod, sdp = list()))

    #niter must be integer
    expect_error(dapper_sample(data_model = dmod,
                               sdp = c(1,2),
                               niter = 3.1))

})

test_that("basic sampler check", {

    posterior_f <- function(dmat, theta) 1
    latent_f    <- function(theta) as.matrix(1)
    statistic_f <- function(xi, sdp, i) 1
    mechanism_f <- function(sdp, sx) 1

    dmod <- new_privacy(posterior_f = posterior_f,
                        latent_f    = latent_f,
                        mechanism_f = mechanism_f,
                        statistic_f = statistic_f,
                        npar = 1)

    out <- dapper_sample(dmod,
                         sdp = 1,
                         init_par = -2,
                         niter = 500)

    expect_equal(out$chain[1], 1)
    expect_equal(out$chain[1], out$chain[250])
})

test_that("return checks work", {

    posterior_f <- function(dmat, theta) 1
    latent_f    <- function(theta) as.matrix(1)
    statistic_f <- function(xi, sdp, i) 1
    mechanism_f <- function(sdp, sx) 1

    #check latent_f()
    latent_e <- function(theta) 1
    dmod <- new_privacy(posterior_f = posterior_f,
                        latent_f    = latent_e,
                        mechanism_f = mechanism_f,
                        statistic_f = statistic_f,
                        npar = 1)

    expect_error(dapper_sample(dmod,
                         sdp = 1,
                         init_par = -2,
                         niter = 500))

    #check posterior_f()
    posterior_e <- function(dmat, theta) as.matrix(1)
    dmod <- new_privacy(posterior_f = posterior_e,
                        latent_f    = latent_f,
                        mechanism_f = mechanism_f,
                        statistic_f = statistic_f,
                        npar = 1)

    expect_error(dapper_sample(dmod,
                               sdp = 1,
                               init_par = -2,
                               niter = 500))

    #check sdp and statistic_f()
    dmod <- new_privacy(posterior_f = posterior_f,
                        latent_f    = latent_f,
                        mechanism_f = mechanism_f,
                        statistic_f = statistic_f,
                        npar = 1)
    
    expect_error(dapper_sample(dmod,
                               sdp = as.matrix(1),
                               init_par = -2,
                               niter = 500))



})

