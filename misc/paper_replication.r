

post_smpl <- function(dmat, theta) {
  x <- cbind(1,dmat[,-1])
  y <- dmat[,1]

  ps_s2 <- solve((1/2) * t(x) %*% x + (1/4) * diag(3))
  ps_m <- ps_s2 %*% (t(x) %*% y) * (1/2)

  MASS::mvrnorm(1, mu = ps_m, Sigma = ps_s2)
}


lik_smpl <- function(theta) {
  #xmat <- MASS::mvrnorm(1 , mu = c(.9,-1.17), Sigma = diag(2))
  xmat <- c(rnorm(1, mean = .9), rnorm(1, mean = -1.17))
  y <- c(1,xmat) %*% theta + rnorm(1, sd = sqrt(2))
  c(y,xmat)
}

clamp_data <- function(dmat) {
  pmin(pmax(dmat,-10),10) / 10
}

tstat <- function(dmat) {
  sdp_mat <- clamp_data(dmat)
  ydp <- sdp_mat[,1, drop = FALSE]
  xdp <- cbind(1,sdp_mat[,-1, drop = FALSE])

  s1 <- t(xdp) %*% ydp
  s2 <- t(ydp) %*% ydp
  s3 <- t(xdp) %*% xdp

  ur_s1 <- c(s1)
  ur_s2 <- c(s2)
  ur_s3 <- s3[upper.tri(s3,diag = TRUE)][-1]
  c(ur_s1,ur_s2,ur_s3)
}

st_update <- function(st, xs, xo) {
  st - tstat(t(xo)) + tstat(t(xs))
}

st_init <- function(dmat) {
  tstat(dmat)
}


gen_priv_zt <- function(epsilon) {
  function(sdp, zt) {
    sum(VGAM::dlaplace(sdp - zt, 0, deltaa/epsilon, TRUE))
  }
}

deltaa <- 13
n <- 5
epsilon <- 5pro
xmat <- MASS::mvrnorm(n, mu = c(.9,-1.17), Sigma = diag(2))
beta <- c(-1.79, -2.89, -0.66)
y <- cbind(1,xmat) %*% beta + rnorm(n, sd = sqrt(2))
z <- tstat(cbind(y,xmat))
z <- z + VGAM::rlaplace(length(z), location = 0, scale = deltaa/epsilon)

dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = gen_priv_zt(epsilon),
                    st_update = st_update,
                    st_calc = st_init,
<<<<<<< HEAD
                    npar = 3)
profvis::profvis({
=======
                    npar = 4)

>>>>>>> 7533517b6202d31970f97642c575591b7e7532ae
tmp <- mcmc_privacy(dmod,
                    sdp = z,
                    nobs = n,
                    init_par = beta,
<<<<<<< HEAD
                    niter = 5000,
=======
                    niter = 10,
>>>>>>> 7533517b6202d31970f97642c575591b7e7532ae
                    chains = 1,
                    varnames = c("beta0", "beta1", "beta2"))})

posterior::summarize_draws(tmp$chain)
bayesplot::mcmc_trace(tmp$chain)



load("C:\\Users\\kevin\\Downloads\\dataaugmentation-mcmc-differentialprivacy-main\\dataaugmentation-mcmc-differentialprivacy-main\\linearregression_sdpfixxy_N100p2noiselevel1iexp1.RData")
z <- c(sdp[[1]], sdp[[2]], sdp[[3]][upper.tri(sdp[[3]],diag = TRUE)][-1])


deltaa <- 13
set.seed(2)
n <- 100
epsilon <- 10
xmat <- MASS::mvrnorm(n, mu = c(.9,-1.17), Sigma = diag(2))
#beta <- runif(3, -4, 4)
beta <- c(-1.79, -2.89, -0.66)

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
                    niter =  3,
                    warmup = 100,
                    chains = 1)


tmp <- sapply(1:10, function(s) lik_smpl(beta)) %>% t()
