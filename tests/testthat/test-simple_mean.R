test_that("test on simple means model", {
    #prior ~ N(0,2ˆ2)
    post_smpl <- function(dmat, theta) {
      x <- c(dmat)
      xbar <- mean(x)
      n <- length(x)
      pr_m <- 0
      pr_s2 <- 4
      ps_s2 <- 1/(1/pr_s2 + n)
      ps_m <- ps_s2 * ((1/pr_s2)*pr_m + n * xbar)
      rnorm(1, mean = ps_m, sd = sqrt(ps_s2))
    }
    clamp_data <- function(dmat) {
      a <- -10
      b <- 10
      if(length(dmat) == 1) {
        min(max(dmat,-10),10)
      } else {
        pmin(pmax(dmat,-10),10)
      }
    }
    tstat <- function(dmat) {
      sum(clamp_data(dmat))
    }
    st_update <- function(st, xs, xo) {
      st - tstat(xo) + tstat(xs)
    }
    st_calc <- function(dmat) {
      tstat(dmat)
    }
    lik_smpl <- function(theta) {
      rnorm(1,mean = theta, sd = 1)
    }
    gen_priv <- function(epsilon) {
      function(sdp, zt) {
        sum(VGAM::dlaplace(sdp - zt, 0, 1/epsilon, TRUE))
      }
    }

    set.seed(1)
    n <- 100
    epsilon <- 20
    y <- rnorm(n, mean = -2)
    y <- sum(clamp_data(y))
    z <- y + extraDistr::rlaplace(1, sigma = 1/epsilon)

    dmod <- new_privacy(post_smpl = post_smpl,
                        lik_smpl = lik_smpl,
                        ll_priv_mech = gen_priv(epsilon),
                        st_update = st_update,
                        st_calc = st_calc,
                        npar = 1)

    tmp <- mcmc_privacy(dmod,
                        sdp = z,
                        nobs = n,
                        init_par = -2,
                        niter = 500)

    expect_equal(summary(tmp$chain)$mean, -2, tolerance = .3)
})

test_that("test single observation", {
  #prior ~ N(0,2ˆ2)
  post_smpl <- function(dmat, theta) {
    x <- c(dmat)
    xbar <- mean(x)
    n <- length(x)
    pr_m <- 0
    pr_s2 <- 4
    ps_s2 <- 1/(1/pr_s2 + n)
    ps_m <- ps_s2 * ((1/pr_s2)*pr_m + n * xbar)
    rnorm(1, mean = ps_m, sd = sqrt(ps_s2))
  }
  clamp_data <- function(dmat) {
    a <- -10
    b <- 10
    if(length(dmat) == 1) {
      min(max(dmat,-10),10)
    } else {
      pmin(pmax(dmat,-10),10)
    }
  }
  tstat <- function(dmat) {
    sum(clamp_data(dmat))
  }
  st_update <- function(st, xs, xo) {
    st - tstat(xo) + tstat(xs)
  }
  st_calc <- function(dmat) {
    tstat(dmat)
  }
  lik_smpl <- function(theta) {
    rnorm(1,mean = theta, sd = 1)
  }
  gen_priv <- function(epsilon) {
    function(sdp, zt) {
      ep <- 1/epsilon
      -abs(sdp - zt)/ep
    }
  }

  set.seed(1)
  n <- 100
  epsilon <- 20
  y <- rnorm(n, mean = -2)
  y <- sum(clamp_data(y))
  z <- y + extraDistr::rlaplace(1, sigma = 1/epsilon)

  dmod <- new_privacy(post_smpl = post_smpl,
                      lik_smpl = lik_smpl,
                      ll_priv_mech = gen_priv(epsilon),
                      st_update = st_update,
                      st_calc = st_calc,
                      npar = 1)

  tmp <- mcmc_privacy(dmod,
                      sdp = z,
                      nobs = n,
                      init_par = -2,
                      niter = 500)

  expect_equal(summary(tmp$chain)$mean, -2, tolerance = .3)
})

test_that("correctly defaults to st_calc with null update", {
    #prior ~ N(0,2ˆ2)
    post_smpl <- function(dmat, theta) {
      x <- c(dmat)
      xbar <- mean(x)
      n <- length(x)
      pr_m <- 0
      pr_s2 <- 4
      ps_s2 <- 1/(1/pr_s2 + n)
      ps_m <- ps_s2 * ((1/pr_s2)*pr_m + n * xbar)
      rnorm(1, mean = ps_m, sd = sqrt(ps_s2))
    }
    clamp_data <- function(dmat) {
      a <- -10
      b <- 10
      if(length(dmat) == 1) {
        min(max(dmat,-10),10)
      } else {
        pmin(pmax(dmat,-10),10)
      }
    }
    tstat <- function(dmat) {
      sum(clamp_data(dmat))
    }

    st_calc <- function(dmat) {
      tstat(dmat)
    }
    lik_smpl <- function(theta) {
      rnorm(1,mean = theta, sd = 1)
    }
    gen_priv <- function(epsilon) {
      function(sdp, zt) {
        ep <- 1/epsilon
        -abs(sdp - zt)/ep
      }
    }

    set.seed(1)
    n <- 100
    epsilon <- 20
    y <- rnorm(n, mean = -2)
    y <- sum(clamp_data(y))
    z <- y + extraDistr::rlaplace(1, sigma = 1/epsilon)

    dmod <- new_privacy(post_smpl = post_smpl,
                        lik_smpl = lik_smpl,
                        ll_priv_mech = gen_priv(epsilon),
                        st_update = NULL,
                        st_calc = st_calc,
                        npar = 1)

    tmp <- mcmc_privacy(dmod,
                        sdp = z,
                        nobs = n,
                        init_par = -2,
                        niter = 500)

    expect_equal(summary(tmp$chain)$mean, -2, tolerance = .3)
})
