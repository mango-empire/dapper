test_that("test on simple means model", {
    #prior ~ N(0,2ˆ2)
    post_f <- function(dmat, theta) {
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
    st_f <- function(dmat) {
      tstat(dmat)
    }
    lik_f <- function(theta) {
      matrix(rnorm(100, mean = theta, sd = 1), ncol = 1)
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
    sdp <- y + extraDistr::rlaplace(1, sigma = 1/epsilon)

    dmod <- new_privacy(post_f = post_f,
                        lik_f = lik_f,
                        priv_f = gen_priv(epsilon),
                        st_f = st_f,
                        add = TRUE,
                        npar = 1)

    tmp <- dapper_sample(dmod,
                        sdp = sdp,
                        init_par = -2,
                        niter = 500)

    expect_equal(summary(tmp$chain)$mean, -2, tolerance = .3)
})

test_that("test single observation", {
  #prior ~ N(0,2ˆ2)
  post_f <- function(dmat, theta) {
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
  st_f <- function(dmat) {
    tstat(dmat)
  }
  lik_f <- function(theta) {
    matrix(rnorm(100, mean = theta, sd = 1), ncol = 1)
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
  sdp <- y + extraDistr::rlaplace(1, sigma = 1/epsilon)

  dmod <- new_privacy(post_f = post_f,
                      lik_f = lik_f,
                      priv_f = gen_priv(epsilon),
                      st_f = st_f,
                      npar = 1)

  tmp <- dapper_sample(dmod,
                    sdp = sdp,
                    init_par = -2,
                    niter = 500)

  expect_equal(summary(tmp$chain)$mean, -2, tolerance = .3)
})

