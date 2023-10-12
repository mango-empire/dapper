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

tstat <- function(dmat) {
  sum(dmat)
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
n <- 16
epsilon <- .1
y <- rnorm(n, mean = -2)
y <- sum(y)
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
                    niter = 15000)

summarise_draws(tmp$chain)

bayesplot::mcmc_trace(tmp$chain)

acf(tmp$chain)

###


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

tstat <- function(dmat) {
  sum(dmat)
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

ll_prop <- function(dmat, theta, theta_prop) {
  dnorm(theta_prop, mean = 0, sd = 2 , log = TRUE)
}

ll_prior <- function(theta) {
  dnorm(theta, mean = 0, sd = 2, log = TRUE)
}

prop_smpl <- function(dmat, theta) {
  rnorm(1, mean = 0, sd = 2)
}


set.seed(1)
n <- 16
epsilon <- .01
y <- rnorm(n, mean = -2)
y <- sum(clamp_data(y))
z <- y + extraDistr::rlaplace(1, sigma = 1/epsilon)

dmod <- new_privacy(prop_smpl = prop_smpl,
                    lik_smpl = lik_smpl,
                    ll_lik = ll_lik,
                    ll_prior = ll_prior,
                    ll_prop = ll_prop,
                    ll_priv_mech = gen_priv(epsilon),
                    st_update = st_update,
                    st_calc = st_calc,
                    npar = 1)

tmp <- mcmc_privacy_prop(dmod,
                         sdp = z,
                         nobs = n,
                         init_par = -2,
                         niter = 5000)

summarise_draws(tmp$chain)


bayesplot::mcmc_trace(tmp$chain)

acf(tmp$chain)

