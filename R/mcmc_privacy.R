mcmc_privacy_add <- function(post_smpl,
                             lik_smpl,
                             ll_priv_mech,
                             st_update,
                             st_init,
                             sdp,
                             npar,
                             nobs,
                             init_par,
                             niter) {
  d_mat <- lapply(1:nobs, function(s) lik_smpl(init_par))
  d_mat <- do.call(rbind, d_mat)
  theta_mat <- matrix(0, nrow = niter, ncol = npar)
  theta <- init_par
  pb <- txtProgressBar(1, niter * n, style=3)
  st <- st_init(d_mat)
  counter <- 0
  for(i in 1:niter) {
    theta_mat[i,] <- theta
    theta <- post_smpl(d_mat, theta)
    for(j in 1:nobs) {
      xs <- lik_smpl(theta)
      xo <- d_mat[j,]
      sn <- st_update(st, xs, xo)
      a <- exp(ll_priv_mech(sdp, sn) - ll_priv_mech(sdp, st))
      if(runif(1) < min(a,1)) {
        counter <- counter + 1
        d_mat[j,] <- xs
        st <- sn
      }
      setTxtProgressBar(pb, i*n + j)
    }
  }
  print(paste(cat("\n","Accept Prob: "), counter/(nobs*niter)))
  theta_mat
}


mcmc_privacy <- function(pobj,
                         sdp,
                         nobs,
                         init_par,
                         niter) {
  d_mat <- lapply(1:nobs, function(s) lik_smpl(init_par))
  d_mat <- do.call(rbind, d_mat)
  theta_mat <- matrix(0, nrow = niter, ncol = npar)
  theta <- init_par
  pb <- txtProgressBar(1, niter * n, style=3)
  st <- st_init(d_mat)
  counter <- 0
  for(i in 1:niter) {
    theta_mat[i,] <- theta
    theta <- post_smpl(d_mat, theta)
    for(j in 1:nobs) {
      xs <- lik_smpl(theta)
      xo <- d_mat[j,]
      sn <- st_update(st, xs, xo)
      a <- exp(ll_priv_mech(sdp, sn) - ll_priv_mech(sdp, st))
      if(runif(1) < min(a,1)) {
        counter <- counter + 1
        d_mat[j,] <- xs
        st <- sn
      }
      setTxtProgressBar(pb, i*n + j)
    }
  }
  print(paste(cat("\n","Accept Prob: "), counter/(nobs*niter)))
  theta_mat
}

