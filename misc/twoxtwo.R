
post_smpl <- function(dmat, theta) {
  x <- c(dmat)
  alpha <- (1 + x) / (4 + sum(x))
  extraDistr::rdirichlet(1, alpha)
}


lik_smpl <- function(theta) {
  rmultinom(1, 300, theta)
}


gen_priv_zt <- function(epsilon) {
  function(sdp, zt) {
    sum(dnorm(sdp - zt, 0, 1/epsilon, log = TRUE))
  }
}


set.seed(2)
n <- 400
epsilon <- .1
theta <- c(.05, .35, .20, .40)
x <- rmultinom(1, n, theta) %>% t()
xmat <- matrix(x, 2, 2, byrow = TRUE)

BioProbability::odds.ratio(xmat, conf.int = TRUE)

z <- z + VGAM::rlaplace(length(z), location = 0, scale = deltaa/epsilon)

dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = gen_priv_zt(epsilon),
                    st_update = st_update,
                    st_calc = st_init,
                    npar = 3)

tmp <- mcmc_privacy(dmod,
                    sdp = z,
                    nobs = n,
                    init_par = beta,
                    niter = 10000,
                    chains = 1,
                    varnames = c("beta0", "beta1", "beta2"))
