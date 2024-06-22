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

    post_f <- function(dmat, theta) 1
    latent_f <- function(theta) as.matrix(1)
    st_f <- function(xi, sdp, i) 1
    priv_f <- function(sdp, tx) 1

    dmod <- new_privacy(post_f = post_f,
                        latent_f = latent_f,
                        priv_f = priv_f,
                        st_f = st_f,
                        npar = 1)

    out <- dapper_sample(dmod,
                         sdp = 1,
                         init_par = -2,
                         niter = 500)

    expect_equal(out$chain[1], 1)
    expect_equal(out$chain[1], out$chain[250])
})
