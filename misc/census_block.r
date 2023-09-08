#determine performance using census block group as upper level invariant
#county = 003
#tract = 030400
#block group = 4
#there are 14 blocks within this block group

m_tot <- c(16, 17, 72, 58, 20, 25, 42, 44, 46, 21, 54, 33, 30, 14)
f_tot <- c(46, 31, 63, 29, 28, 83, 69, 42, 70, 49, 35, 15, 65, 34)

post_smpl <- function(dmat, theta) {
  S <- c(dmat)
  bo <- S[as.logical((1:28) %% 2)]
  be <- S[!as.logical((1:28) %% 2)]
  stot <- bo + be
  N <- sum(stot)
  alpha <- dirmult::rdirichlet(1, stot + 1)
  theta1 <- rbeta(length(bo), bo + 1, be + 1)
  c(alpha, theta1)
}

lik_smpl <- function(theta) {
  N <- 1151
  alpha <- theta[1:14]
  theta1 <- theta[-c(1:14)]
  alpha <- dirmult::rdirichlet(1, rep(1,14))
  stot <- c(rmultinom(1, N, alpha))
  theta <- runif(14)
  b1 <- rbinom(14, stot, theta)
  b2 <- stot - b1

  c(rbind(b1,b2))
}


gen_priv <- function(epsilon) {
  function(sdp, zt) {
    sum(dnorm(sdp-zt, mean = 0, sd = 1/epsilon, log = TRUE))
  }
}

st_calc <- function(dmat) {
  dmat
}

#test posterior sample
set.seed(1)
alpha <- dirmult::rdirichlet(1, rep(1,14))
stot <- c(rmultinom(1, 3000000, alpha))
theta <- runif(14)
b1 <- rbinom(14, stot, theta)
b2 <- stot - b1

dmat <- c(rbind(b1,b2))


tmp <- post_smpl(dmat,1)
tmp[1:14]
tmp[-c(1:14)]



B <- c(rbind(m_tot,f_tot))
sigma <- 26.64 #block group
#sigma <- 200
dp_noise <-  rnorm(28, mean = 0, sd = sigma)
dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = gen_priv(1/sigma),
                    st_update = NULL,
                    st_calc = st_calc,
                    npar = 28)

vn <- c(paste0('alpha',1:14), paste0('theta',1:14))

tmp <- mcmc_privacy(dmod,
                    sdp = B + dp_noise,
                    nobs = 1,
                    init_par = rep(.5, 28),
                    niter = 10000,
                    chains = 1,
                    varnames = vn)

summary.dpsim(tmp)

stab <- summary(tmp)
stab

abs(stab$q5 - stab$q95) * 100


tdf <- sapply(1:10000, function(s) post_smpl(B, 1)) %>% t() %>% posterior::as_draws_matrix()
tsf <- summary(tdf)

abs(tsf$q5 - tsf$q95) * 100


#independent metropolis

m_tot <- c(16, 17, 72, 58, 20, 25, 42, 44, 46, 21, 54, 33, 30, 14)
f_tot <- c(46, 31, 63, 29, 28, 83, 69, 42, 70, 49, 35, 15, 65, 34)

post_smpl <- function(dmat, theta) {
  S <- c(dmat)
  bo <- S[as.logical((1:28) %% 2)]
  be <- S[!as.logical((1:28) %% 2)]
  stot <- bo + be
  N <- sum(stot)
  alpha <- dirmult::rdirichlet(1, stot + 1)
  theta1 <- rbeta(length(bo), bo + 1, be + 1)
  c(alpha, theta1)
}

lik_smpl <- function(theta) {
  alpha <- dirmult::rdirichlet(1, rep(1,14))
  stot <- c(rmultinom(1, 3000000, alpha))
  theta <- runif(14)
  b1 <- rbinom(14, stot, theta)
  b2 <- stot - b1
  c(rbind(b1,b2))
}

lik_smpl <- function(theta) {
  c(rbind(m_tot,f_tot))
}


gen_priv <- function(epsilon) {
  function(sdp, zt) {
    sum(dnorm(sdp-zt, mean = 0, sd = 1/epsilon, log = TRUE))
  }
}

st_calc <- function(dmat) {
  dmat
}

#test posterior sample
set.seed(1)
alpha <- dirmult::rdirichlet(1, rep(1,14))
stot <- c(rmultinom(1, 3000000, alpha))
theta <- runif(14)
b1 <- rbinom(14, stot, theta)
b2 <- stot - b1

dmat <- c(rbind(b1,b2))


tmp <- post_smpl(dmat,1)
tmp[1:14]
tmp[-c(1:14)]



B <- c(rbind(m_tot,f_tot))
sigma <- 26.64 #block group
#sigma <- 200
dp_noise <-  rnorm(28, mean = 0, sd = sigma)
dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = gen_priv(1/sigma),
                    st_update = NULL,
                    st_calc = st_calc,
                    npar = 28)

vn <- c(paste0('alpha',1:14), paste0('theta',1:14))

tmp <- mcmc_privacy(dmod,
                    sdp = B + dp_noise,
                    nobs = 1,
                    init_par = rep(.5, 28),
                    niter = 50000,
                    chains = 1,
                    varnames = vn)

stab <- summary(tmp)
stab


#independent metrop

m_tot <- c(16, 17, 72, 58, 20, 25, 42, 44, 46, 21, 54, 33, 30, 14)
f_tot <- c(46, 31, 63, 29, 28, 83, 69, 42, 70, 49, 35, 15, 65, 34)

post_smpl <- function(dmat, theta) {
  S <- c(dmat)
  bo <- S[as.logical((1:28) %% 2)]
  be <- S[!as.logical((1:28) %% 2)]
  stot <- bo + be
  N <- sum(stot)
  alpha <- dirmult::rdirichlet(1, stot + 1)
  theta1 <- rbeta(length(bo), bo + 1, be + 1)
  c(alpha, theta1)
}

lik_smpl <- function(theta) {
  N <- 1151
  alpha <- dirmult::rdirichlet(1, rep(1,14))
  stot <- c(rmultinom(1, N, alpha))
  theta <- runif(14)
  b1 <- rbinom(14, stot, theta)
  b2 <- stot - b1

  dmat <- c(rbind(b1,b2))
}


gen_priv <- function(epsilon) {
  function(sdp, zt) {
    sum(dnorm(sdp-zt, mean = 0, sd = 1/epsilon, log = TRUE))
  }
}

st_calc <- function(dmat) {
  dmat
}



B <- c(rbind(m_tot,f_tot))
sigma <- 26.64 #block group
#sigma <- 200
dp_noise <-  rnorm(28, mean = 0, sd = sigma)
dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = gen_priv(1/sigma),
                    st_update = NULL,
                    st_calc = st_calc,
                    npar = 28)

vn <- c(paste0('alpha',1:14), paste0('theta',1:14))

tmp <- mcmc_privacy(dmod,
                    sdp = B + dp_noise,
                    nobs = 1,
                    init_par = rep(.5, 28),
                    niter = 50000,
                    chains = 1,
                    varnames = vn)
