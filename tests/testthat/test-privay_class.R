test_that("inputs are verified", {
  post_f <- function(dmat, theta) 1
  latent_f <- function(theta) as.matrix(1)
  st_f <- function(xi, sdp, i) 1
  priv_f <- function(sdp, tx) 1

  #names
  post_e <- function(mat, theta) NULL
  expect_error(new_privacy(post_f = post_e,
                           latent_f = latent_f,
                           priv_f = priv_f,
                           st_f = st_f,
                           npar = 1))

  latent_e <- function(value) NULL
  expect_error(new_privacy(post_f = post_f,
                           latent_f = latent_e,
                           priv_f = priv_f,
                           st_f = st_f,
                           npar = 1))

  st_e <- function(xi, sdp, j) NULL
  expect_error(new_privacy(post_f = post_f,
                           latent_f = latent_f,
                           priv_f = priv_f,
                           st_f = st_e,
                           npar = 1))

  ste <- function(xi, sdp) NULL
  expect_error(new_privacy(post_f = post_f,
                           latent_f = latent_f,
                           priv_f = priv_f,
                           st_f = st_e,
                           npar = 1))

  priv_e <- function(sdp, te) NULL
  expect_error(new_privacy(post_f = post_f,
                           latent_f = latent_f,
                           priv_f = priv_e,
                           st_f = st_f,
                           npar = 1))



})
