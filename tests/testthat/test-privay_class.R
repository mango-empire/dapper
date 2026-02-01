test_that("inputs are verified", {
  post_f   <- function(dmat, theta) 1
  latent_f <- function(theta) as.matrix(1)
  stat_f   <- function(xi, sdp, i) 1
  mech_f   <- function(sdp, tx) 1

  #names
  post_e <- function(mat, theta) NULL
  expect_error(new_privacy(posterior_f = post_e,
                           latent_f    = latent_f,
                           mechanism_f = mech_f,
                           statistic_f = stat_f,
                           npar = 1))

  latent_e <- function(value) NULL
  expect_error(new_privacy(posterior_f = post_f,
                           latent_f    = latent_e,
                           mechanism_f = mech_f,
                           statistic_f = stat_f,
                           npar = 1))

  stat_e <- function(xi, sdp, j) NULL
  expect_error(new_privacy(posterior_f = post_f,
                           latent_f    = latent_f,
                           mechanism_f = mech_f,
                           statistic_f = stat_e,
                           npar = 1))

  mech_e <- function(sdp, te) NULL
  expect_error(new_privacy(posterior_f = post_f,
                           latent_f    = latent_f,
                           mechanism_f = mech_e,
                           statistic_f = stat_f,
                           npar = 1))



})
