library(mcmc)
library(adaptMCMC)

lupost_factory <- function(x, y) function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  mux <- theta[3]
  sx <- exp(theta[4])
  sy <- exp(theta[5])

  t1 <- dnorm(x, mean = mux, sd = exp(sx), log = TRUE)
  t2 <- dnorm(y, mean = alpha + beta * x, sd = exp(sy), log = TRUE)

  sum(t1) + sum(t2)
}

set.seed(1)
n <- 100
alpha <- 2
beta <- -3
x <- rnorm(n, 1, 3)
y <- alpha + beta*x + rnorm(n,0,5)
#----------

lupost_factory <- function(x, y) function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  sy <- exp(theta[3])
  t1 <- dnorm(y, mean = alpha + beta * x, sd = sy, log = TRUE)
  sum(t1)
}

post_smpl <- function(dmat, theta) {
  x <- dmat[,1]
  y <- dmat[,2]
  mout <- metrop(lupost_factory(x,y),
                 initial = theta,
                 scale = c(.5,.5,.1),
                 nbatch = 1)
  c(mout$batch)
}

post_smpl <- function(dmat, theta) {
  x <- dmat[,1]
  y <- dmat[,2]
  invisible(capture.output(
  mout <- adaptMCMC::MCMC(lupost_factory(x,y),
                          n= 200,
                          init= theta,
                          scale=c(1, 1 ,.1),
                          adapt=TRUE,
                          list = FALSE,
                          showProgressBar = FALSE,
                          acc.rate=0.234)))
  c(mout[200,])
}



lik_smpl <- function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  sy <- exp(theta[3])
  x <- runif(1)
  y <- rnorm(1, mean = alpha + beta*x, sd = sy)
  c(x,y)
}

st_calc <- function(dmat) {
  x <- dmat[,1]
  y <- dmat[,2]
  beta <- cov(x,y)/var(x)
  alpha <- mean(y) - beta * mean(x)
  n <- length(y)
  c((n-1) * cov(x,y), (n-1) * var(x), alpha)
}

priv_mech_factory <- function(epsilon, delta12, delta3) {
  function(sdp, xt) {
    #sum(VGAM::dlaplace(sdp - xt, 0, 3/epsilon, TRUE))
    t1 <- VGAM::dlaplace(sdp[1] - xt[1], 0, 3*delta12/epsilon, TRUE)
    t2 <- VGAM::dlaplace(sdp[2] - xt[2], 0, 3*delta12/epsilon, TRUE)
    t3 <- VGAM::dlaplace(sdp[3] - xt[3], 0, 3*delta3/epsilon, TRUE)
    sum(t1 + t2 + t3)
  }
}


set.seed(1)
epsilon <- 3
n <- 50
delta12 <- 1 - 1/n
alpha <- 5
beta <- -3
x <- runif(n)
y <- alpha + beta*x + rnorm(n,0,3)
sdp <- st_calc(cbind(x,y))



sdp[1] <- sdp[1] + VGAM::rlaplace(1, 0, delta12 * 3/epsilon)
sdp[2] <- sdp[2] + VGAM::rlaplace(1, 0, delta12 * 3/epsilon)
beta_sdp <- sdp[1]/sdp[2]
delta3 <- (1/n) * (1 + abs(beta_sdp))
sdp[3] <- mean(y) - beta_sdp * mean(x) + VGAM::rlaplace(1, 0, delta3 * 3/epsilon)


dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = priv_mech_factory(epsilon,delta12,delta3),
                    st_calc = st_calc,
                    add = FALSE,
                    npar = 3)

tmp <- mcmc_privacy(dmod,
               sdp = sdp,
               nobs = n,
               init_par = c(1,1,1),
               niter = 10000,
               chains = 1,
               varnames = c("alpha", "beta", "sigma"))

summary(tmp)
bayesplot::mcmc_trace(tmp$chain)
 #test posterior sampler

p.log <- lupost_factory(x,y)

samp <- MCMC(p.log, n=2000, init=c(1, 1, 1), scale=c(1, 1 ,1),
             adapt=TRUE, acc.rate=0.234, list=FALSE)


#-------------------------------------

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

standata <- list(N=n, x = x, y= y)


stanmodel1 <- stan_model("misc/alabi/linear.stan")

reg_mod <- cmdstan_model("misc/alabi/linear.stan")
reg_fit <- reg_mod$sample(data = standata, chains = 1)

ed <- diag(reg_fit$inv_metric()$'1')

reg_fit2 <- reg_mod$sample(data = standata,
                           chains = 1,
                           adapt_engaged = FALSE,
                           metric = "diag_e",
                           inv_metric = ed,
                           step_size = .3,
                           adapt_delta = .6)



reg_fit2$draws()
samples1 <- sampling(stanmodel1, data = standata, init = c(1,1,1),
                     iter = 15, warmup = 10, chains = 1)
extract(samples1)

samples1 <- sampling(stanmodel1, data = standata, init = c(1,1,1),
                     iter = 5000, warmup = 1000, chains = 1,
                     verbose=FALSE, refresh =0)

 tmp <- get_sampler_params(samples1)[[1]]


post_smpl <- function(dmat, theta) {
  x <- dmat[,1]
  y <- dmat[,2]
  standata <- list(N=n, x = x, y= y)
  suppressWarnings(
  samples1 <- sampling(stanmodel1, data = standata, init = theta,
                       iter = 100, warmup = 50, chains = 1, refresh=0))
  tmp <- extract(samples1)
  alpha <- tmp$alpha[50]
  beta <- tmp$beta[50]
  sigma <- tmp$sigma[50]
  c(alpha,beta, sigma)
}

lik_smpl <- function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  sy <- theta[3]
  x <- runif(1)
  y <- rnorm(1, mean = alpha + beta*x, sd = sy)
  c(x,y)
}


dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = priv_mech_factory(epsilon),
                    st_calc = st_calc,
                    add = FALSE,
                    npar = 3)

tmp <- mcmc_privacy(dmod,
             sdp = sdp,
             nobs = n,
             init_par = c(1,1,1),
             niter = 1500,
             warmup = 500,
             chains = 1,
             varnames = c("alpha", "beta", "sigma"))


#-------------------------------------------------------------------------------


lik_smpl <- function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  x <- runif(1)
  y <- rnorm(1, mean = alpha + beta * x, sd = 3)
  c(x,y)
}

post_smpl <- function(dmat, theta) {
  x <- dmat[,1]
  y <- dmat[,2]
  xm <- cbind(1, x)
  Si <- (1/9) * t(xm) %*% xm + (1/100) * diag(2)
  mu <- (1/9) * solve(Si) %*% t(xm) %*% y
  MASS::mvrnorm(1, mu = mu, Sigma = solve(Si))
}

st_calc <- function(dmat) {
  x <- dmat[,1]
  y <- dmat[,2]
  n <- length(y) - cov(x,y)/var(x)
  s1 <- (n-1) * cov(x,y)
  s2 <- (n-1) * var(x)
  s3 <- mean(y)
  s4 <- mean(x)
  c(s1, s2, s3, s4)
}

priv_mech_factory <- function(n, epsilon) {
  function(sdp, xt) {
    delta1 <- (1- 1/n)
    delta3 <- 1/n
    t1 <- VGAM::dlaplace(sdp[1] - xt[1], 0, 3 * delta1/epsilon, TRUE) #EXPENSIVE
    t2 <- VGAM::dlaplace(sdp[2] - xt[2], 0, 3 * delta1/epsilon, TRUE)
    t3 <- VGAM::dlaplace(sdp[3] - xt[3], 0, 3 * delta3/epsilon, TRUE)
    t4 <- VGAM::dlaplace(sdp[4] - xt[4], 0, 3 * delta3/epsilon, TRUE)
    sum(c(t1,t2,t3,t4))
  }
}



set.seed(1)
epsilon <- 1
n <- 50
alpha <- 5
beta <- -3
x <- runif(n)
y <- alpha + beta*x + rnorm(n,0,3)
delta1 <- 1 - 1/n
delta3 <- 1/n
sdp <- st_calc(cbind(x,y))
sdp[1] <- sdp[1] + VGAM::rlaplace(1, 0, 3 * delta1/epsilon)
sdp[2] <- sdp[2] + VGAM::rlaplace(1, 0, 3 * delta1/epsilon)
sdp[3] <- sdp[3] + VGAM::rlaplace(1, 0, 3 * delta3/epsilon)
sdp[4] <- sdp[4] + VGAM::rlaplace(1, 0, 3 * delta3/epsilon)



dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = priv_mech_factory(n, epsilon),
                    st_calc = st_calc,
                    add = FALSE,
                    npar = 2)
handlers("cli")
with_progress({
tmp <- gdp_sample(dmod,
                  sdp = sdp,
                  nobs = n,
                  init_par = c(1,1),
                  niter = 500,
                  varnames = c("alpha", "beta"))})

summary(tmp)


#-------------------------------------------------------------------------------



lik_smpl <- function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  sy <- exp(theta[3])
  x <- runif(1)
  y <- rnorm(1, mean = alpha + beta * x, sd = sy)
  c(x,y)
}

lupost_factory <- function(x, y) {
  function(theta) {
    alpha <- theta[1]
    beta <- theta[2]
    sy <- exp(theta[3])
    t1 <- dnorm(theta[3], mean =  0, sd = 10, log = TRUE)
    t2 <- dnorm(y, mean = alpha + beta * x, sd = sy, log = TRUE)
    sum(c(t1,t2))
  }
}

post_smpl <- function(dmat, theta) {
  x <- dmat[,1]
  y <- dmat[,2]
  invisible(capture.output(
    mout <- adaptMCMC::MCMC(lupost_factory(x,y),
                            n= 200,
                            init= theta,
                            scale=c(1, 1 ,.1),
                            adapt=TRUE,
                            list = FALSE,
                            showProgressBar = FALSE,
                            acc.rate=0.234)))
  c(mout[200,])
}

st_calc <- function(dmat) {
  x <- dmat[,1]
  y <- dmat[,2]
  n <- length(y) - cov(x,y)/var(x)
  s1 <- (n-1) * cov(x,y)
  s2 <- (n-1) * var(x)
  s3 <- mean(y)
  s4 <- mean(x)
  c(s1, s2, s3, s4)
}

priv_mech_factory <- function(n, epsilon) {
  function(sdp, xt) {
    delta1 <- (1- 1/n)
    delta3 <- 1/n
    t1 <- VGAM::dlaplace(sdp[1] - xt[1], 0, 3 * delta1/epsilon, TRUE) #EXPENSIVE
    t2 <- VGAM::dlaplace(sdp[2] - xt[2], 0, 3 * delta1/epsilon, TRUE)
    t3 <- VGAM::dlaplace(sdp[3] - xt[3], 0, 3 * delta3/epsilon, TRUE)
    t4 <- VGAM::dlaplace(sdp[4] - xt[4], 0, 3 * delta3/epsilon, TRUE)
    sum(c(t1,t2,t3,t4))
  }
}



set.seed(1)
epsilon <- 1
n <- 50
alpha <- 5
beta <- -3
x <- runif(n)
y <- alpha + beta*x + rnorm(n,0,3)
delta1 <- 1 - 1/n
delta3 <- 1/n
sdp <- st_calc(cbind(x,y))
sdp[1] <- sdp[1] + VGAM::rlaplace(1, 0, 3 * delta1/epsilon)
sdp[2] <- sdp[2] + VGAM::rlaplace(1, 0, 3 * delta1/epsilon)
sdp[3] <- sdp[3] + VGAM::rlaplace(1, 0, 3 * delta3/epsilon)
sdp[4] <- sdp[4] + VGAM::rlaplace(1, 0, 3 * delta3/epsilon)



dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = priv_mech_factory(n, epsilon),
                    st_calc = st_calc,
                    add = FALSE,
                    npar = 3)

tmp <- mcmc_privacy(dmod,
                      sdp = sdp,
                      nobs = n,
                      init_par = c(1,1,1),
                      niter = 5000,
                      chains = 1,
                      varnames = c("alpha", "beta", "logsigma"))

summary(tmp)

bayesplot::mcmc_trace(tmp$chain)





