#prior ~ N(0,2ˆ2)
post_smpl <- function(dmat, theta) {
  x <- c(dmat)
  xbar <- mean(x)
  n <- length(x)
  pr_m <- 0
  pr_s2 <- 100
  ps_s2 <- 1/(1/pr_s2 + n)
  ps_m <- ps_s2 * ((1/pr_s2)*pr_m + n * xbar)
  rnorm(1, mean = ps_m, sd = sqrt(ps_s2))
}
tstat <- function(dmat) {
  dmat
}
st_calc <- function(dmat) {
  tstat(dmat)
}
lik_smpl <- function(theta) {
  rnorm(1,mean = theta, sd = 1)
}
gen_priv <- function(epsilon) {
  function(sdp, zt) {
    dnorm(sdp, mean = zt, sd = 1/epsilon, log = TRUE)
  }
}

set.seed(1)
n <- 1
epsilon <- .001
y <- rnorm(n, mean = -2)
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
                    init_par = 40,
                    niter = 1000)

bayesplot::mcmc_trace(tmp$chain)
posterior::summarise_draws(tmp$chain)
