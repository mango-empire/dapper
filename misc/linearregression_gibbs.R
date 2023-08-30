
rm(list = ls());
N <-  100;
p <- 2;
iexp <- 1;
noiselevel <- 1;
niter <- 10000; ## number of mcmc iterations


### add the Rscript read commands here
invisible(eval(parse(text=commandArgs(TRUE))));
set.seed(iexp); ## random seed

## load data generating parameters-------
load( paste("misc/linearregression_datageneratingparam_N", N, "p", p, ".RData", sep = ""));
A <- chol(Sigma_star);
## load sdp----------
load(paste("misc/linearregression_sdpfixxy_N",N,"p",p,"noiselevel",noiselevel,"iexp", iexp, ".RData", sep = ""));

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

sdp_convert <- function(z) {
  s1 <- t(t(z[1:3]))
  s2 <- z[4]
  s3 <- matrix(0, ncol = 3, nrow = 3)
  s3[upper.tri(s3,diag = TRUE)] <- c(100, z[5:9])
  s3[lower.tri(s3)] <- t(s3)[lower.tri(s3)]
  list(s1,s2,s3)
}


set.seed(2)
n <- 100
epsilon <- 0.1
eps <- 0.1
xmat <- MASS::mvrnorm(n, mu = c(.9,-1.17), Sigma = diag(2))
beta <- c(-1.79, -2.89, -0.66)
y <- cbind(1,xmat) %*% beta + rnorm(n, sd = sqrt(2))
z <- tstat(cbind(y,xmat))
z <- z + VGAM::rlaplace(length(z), location = 0, scale = deltaa/epsilon)
sdp <- sdp_convert(z)





## define some helper functions------

## clamping + location-scale transformation
## use this function to normalize each covariate x[-,j] and y
normalize <- function(vec, left, right){
  return(2 * (pmin(right, pmax(left, vec)) - left) / (right - left) - 1 );
}

## initialize the Gibbs sampler----
rinit_gibbs <- function(){
  ## generate beta from prior
  beta <- rnorm(p + 1, 0, tau_star);
  ### generate X from model
  z <- matrix(rnorm(N * p), nrow = p, ncol = N);
  x <- apply(z, 2, function(z_) m_star + A %*% z_); ## Sigma = A A'
  X <- cbind(1, t(x));
  ### generate Y given Xbeta + noise
  Y <- X %*% beta + rnorm(N, 0, sqrt(sigma2_star));
  ### clamp and get diffT1, diffT2, diffT3
  return(state = list(X = X, Y = Y, beta = beta));
}

## compute distance from summary statistics to sdp--------
get_diffT <- function(state, clbounds, sdp){
  ### clamp X and Y
  X_cl <- apply(state$X, 2, function(x) normalize(x, clbounds[1], clbounds[3]));
  X_cl[, 1] <- 1;
  Y_cl <- apply(state$Y, 2, function(y) normalize(y, clbounds[2], clbounds[4]));
  ### compute summary statistics and distance to noisy sdp
  T1 <- t(X_cl) %*% Y_cl - sdp[[1]];
  T2 <- t(Y_cl) %*% Y_cl - sdp[[2]];
  T3 <- t(X_cl) %*% X_cl - sdp[[3]];
  ### save the T(x,y) - sdp
  state$T1 <- T1;
  state$T2 <- T2;
  state$T3 <- T3;
  return(state);
}

## compute density of privacy noise---------
## using Laplace densities
## operations on log-scale for stability
get_logeta <- function(state, sdp){
  return(sum(VGAM::dlaplace(state$T1 , 0, deltaa/ eps, TRUE)) +
    sum(VGAM::dlaplace(state$T2, 0, deltaa/eps, TRUE)) +
      sum(VGAM::dlaplace(state$T3[lower.tri(state$T3)], 0, deltaa/eps, TRUE)) +
      sum(VGAM::dlaplace(diag(state$T3), 0, deltaa / eps, TRUE)));
}


## difference in each individual contribution to t(x,y)
## between current state and proposed state
## ti(xi_star,yi_star) - to(xi, yi)
get_diffeacht <- function(x, x_, y, y_, clbounds){
  xcl <- normalize(x, clbounds[1], clbounds[3]);
  xcl[1] <- 1;
  ycl <- normalize(y, clbounds[2], clbounds[4]);
  xcl_ <- normalize(x_, clbounds[1], clbounds[3]);
  xcl_[1] <- 1;
  ycl_ <- normalize(y_, clbounds[2], clbounds[4]);
  return(diffti = list(diffti1 = xcl_ %o% ycl_ - xcl %o% ycl,
                       diffti2 = ycl_**2 - ycl**2,
                       diffti3 = xcl_ %o% xcl_ - xcl %o% xcl));
}

## independent-MH kernel for xi, yi given xnoti, ynoti, beta, sdp, and other parameters
## propose (xi_star,yi_star) pair from the model f(x,y | theta);
## acceptance ratio is eta(tstar) / eta(t);
update_XiYi_gibbs <- function(i, state, eps){
  acc <- 0;
  ## propose xi and yi from model
  zi <- rnorm(p);
  xi <- c(1, m_star + t(A) %*% zi);
  yi <- xi %*% state$beta + sqrt(sigma2_star) * rnorm(1);
  ## compute T*
  diffti <- get_diffeacht(state$X[i,],xi, state$Y[i], yi, clbounds);
  diffti1_ <- diffti$diffti1; ## xty
  diffti2_ <- diffti$diffti2; ##yty
  diffti3_ <- diffti$diffti3; ## xtx
  logeta_ <- sum(VGAM::dlaplace(state$T1 + diffti1_ , 0, deltaa/ eps, TRUE)) +
    sum(VGAM::dlaplace(state$T2 + diffti2_, 0, deltaa/eps, TRUE)) +
    sum(VGAM::dlaplace((state$T3 + diffti3_)[lower.tri(diffti3_)], 0, deltaa/eps, TRUE)) +
    sum(VGAM::dlaplace(diag(state$T3 + diffti3_), 0, deltaa / eps, TRUE));
  ## accept and reject step
  logu <- log(runif(1));
  if(logu < logeta_ - state$logeta){
    acc <- 1;
    ## accept and change all associated values of x and y
    state$X[i,] <- xi;
    state$Y[i] <- yi;
    state$T1 <- state$T1 + diffti1_;
    state$T2 <- state$T2 + diffti2_;
    state$T3 <- state$T3 + diffti3_;
    state$logeta <- logeta_;
  }
  state$acc <- acc;
  return(state)
}


## update beta given (x,y) using conjugate posterior
update_beta_gibbs <- function(state){
  ## precision: Sigma0^{-1} + 1/sigma^2 * xtx
  ## mean = 1 / precision * (Sigma0^{-1}mu0 + 1 / sigma2 * xty)
  postprec <- t(state$X) %*% state$X / sigma2_star + diag(p+1) / tau_star**2;
  postcov <- solve(postprec);
  postmean <- postcov %*% t(state$X) %*% state$Y / sigma2_star;
  beta <- MASS::mvrnorm(n = 1, postmean, postcov);
  state$beta <- beta;
  return(state);
}


## algorithm 1 in the paper for linear regression
onestep_gibbs <- function(state){
  state <- update_beta_gibbs(state);### step1 of gibbs sampler
  ### step2
  acc <- 0;
  for(i in 1 : N){
    state <- update_XiYi_gibbs(i, state, eps);
    acc <- state$acc + acc;
  }
  state$accMean <- acc / N;
  return(state);
}

diagn <- diag(N);


## initiate the chain
state <- rinit_gibbs();
state <- get_diffT(state, clbounds, sdp);
state$logeta <- get_logeta(state, sdp);

## store parameter values, acceptance rate, and log-posterior
logeta_chain <- double(niter);
beta_chain <- matrix(NA, nrow = p + 1, ncol = niter);
acc_chain <- double(niter);

for(iter in 1 : niter){
  state <- onestep_gibbs(state);
  logeta_chain[iter] <- state$logeta;
  beta_chain[,iter] <- state$beta;
  acc_chain[iter] <- state$accMean;
  if(iter %% 100 == 0) cat("iteration", iter, "\n");
}

plot(logeta_chain[(niter - 2000) : niter], type = "l");
matplot(t(beta_chain[, niter - 2000 : niter ]), type = "l");
summary(acc_chain, type = "l")
