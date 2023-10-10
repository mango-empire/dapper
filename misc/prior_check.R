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
n <- 16
epsilon <- .001
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
                    niter = 15000)

summarise_draws(tmp$chain)

bayesplot::mcmc_trace(tmp$chain)

acf(tmp$chain)

##

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

ll_lik <- function(dmat, theta) {
  x <- c(dmat)
  n <- length(x)
  sum(dnorm(x, mean = theta, sd = 1, log = TRUE))
}

ll_prop <- function(theta, theta_prop) {
  dnorm(theta_prop, mean = 0, sd = 2, log = TRUE)
}

ll_prior <- function(theta) {
  dnorm(theta, mean = 0, sd = 2, log = TRUE)
}

prop_smpl <- function(theta) {
  rnorm(1, mean = 0, sd = 2)
}


set.seed(1)
n <- 16
epsilon <- .001
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
                    niter = 15000)

summarise_draws(tmp$chain)


bayesplot::mcmc_trace(tmp$chain)

acf(tmp$chain)


##fat tailed

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

ll_lik <- function(dmat, theta) {
  x <- c(dmat)
  n <- length(x)
  sum(dnorm(x, mean = theta, sd = 1, log = TRUE))
}

ll_prop <- function(theta, theta_prop) {
  dt(theta_prop, df = 2, log = TRUE)
}

ll_prior <- function(theta) {
  dnorm(theta, mean = 0, sd = 2, log = TRUE)
}

prop_smpl <- function(theta) {
  rt(1, df = 2)
}


set.seed(1)
n <- 16
epsilon <- .001
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
                         niter = 15000)

summarise_draws(tmp$chain)


bayesplot::mcmc_trace(tmp$chain)

acf(tmp$chain)
