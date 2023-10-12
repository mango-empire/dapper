
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

ddp <- function(s, x) {
  xt <- tstat(x)
  exp(sum(VGAM::dlaplace(s - xt, 0, deltaa/epsilon, TRUE)))
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


set.seed(1)
deltaa <- 13
epsilon <- 2
n <- 100
xmat <- MASS::mvrnorm(n, mu = c(.9,-1.17), Sigma = diag(2))
beta <- c(-1.79, -2.89, -0.66)
y <- cbind(1,xmat) %*% beta + rnorm(n, sd = sqrt(2))
sdp <- tstat(cbind(y,xmat))
sdp <- sdp + VGAM::rlaplace(length(sdp), location = 0, scale = deltaa/epsilon)


tmpdf <- as_tibble(cbind(y, xmat))
lm(V1 ~ ., data =tmpdf)

nsim <- 10000
xm <- cbind(y,xmat)
xm <- xm - xm
beta_mat <- matrix(0, nrow = nsim, ncol = 3)


for(j in 1:n) xm[j,] <- lik_smpl(beta_mat[1,])
imp_o <- xm
fh_o <- ddp(sdp, imp_o)

for(i in 2:nsim) {
  beta <- rnorm(3, 0, 2)

  for(j in 1:n) xm[j,] <- lik_smpl(beta)
  fh_n <- ddp(sdp, xm)


  if(log(runif(1)) <= (log(fh_n) - log(fh_o))) {
    beta_mat[i,] <- beta
    fh_o <- fh_n
  } else {
    beta_mat[i,] <- beta_mat[i-1,]
  }
}


tmp <- as_draws(beta_mat[-c(1:(nsim/2)),])
bayesplot::mcmc_trace(tmp)
summarise_draws(tmp)


