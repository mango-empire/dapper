#' Sample from posterior given private data
#'
#' @param data_model A data model created using the privacy class.
#' @param sdp The observed privatized data.
#' @param nobs The number of observations.
#' @param init_par Initial starting point of the chain.
#' @param niter Number of draws.
#' @param warmup Number of iterations to discard as warmup. Default is half of niter.
#'
#' @return A matrix of posterior samples.
#' @export
#'
#' @examples
gdp_sample <- function(data_model,
                       sdp,
                       nobs,
                       init_par,
                       niter = 2000,
                       warmup = floor(niter / 2),
                       chains = 1) {
  #check inputs
  checkmate::qassert(chains, "X?(0,)")
  if(length(init_par) != data_model$npar) stop("Dimension of initial parameter does not match privacy model")

  pb_size <- floor(niter / 100)
  p <- progressr::progressor(pb_size * chains)
  fout <- furrr::future_map(rep(niter, chains), gdp_chain,
                            data_model = data_model,
                            sdp = sdp,
                            init_par = init_par,
                            warmup = warmup,
                            prg_bar = p,
                            .options = furrr::furrr_options(seed = TRUE))

  theta_clist <- lapply(1:chains, function(s) fout[[s]]$sample)
  accept_mat <- do.call(cbind, lapply(1:chains, function(s) fout[[s]]$accept_rate))
  new_gdpout(theta_clist, accept_mat, data_model$varnames)
}




#' single chain sample
#' @export
gdp_chain <- function(data_model,
                      sdp,
                      init_par,
                      niter = 2000,
                      warmup = floor(niter / 2),
                      prg_bar = NULL) {
  #check inputs
  checkmate::qassert(niter, "X?(0,)")
  checkmate::assert_class(data_model, "privacy")

  post_smpl    <- data_model$post_smpl
  lik_smpl     <- data_model$lik_smpl
  ll_priv_mech <- data_model$ll_priv_mech
  st_calc      <- data_model$st_calc
  npar         <- data_model$npar

  accept_rate <- numeric(niter)
  theta_clist <- list()
  dmat        <- lik_smpl(init_par)
  theta_mat   <- matrix(0, nrow = niter, ncol = npar)
  theta       <- init_par
  st          <- st_calc(dmat)
  nobs        <- nrow(dmat)

  for (i in 1:niter) {
    counter        <- 0
    theta_mat[i, ] <- theta
    theta          <- post_smpl(dmat, theta)
    smat           <- lik_smpl(theta)
    for (j in 1:nobs) {
      xs <- smat[j, ]
      xo <- dmat[j, ]
      sn <- NULL
      if (!data_model$add) {
        nmat      <- dmat
        nmat[j, ] <- xs
        sn        <- st_calc(nmat)
      } else {
        #convert xo,xs back to matrices
        sn <- st - st_calc(t(xo)) + st_calc(t(xs))
      }
      a <- exp(ll_priv_mech(sdp, sn) - ll_priv_mech(sdp, st))
      if (stats::runif(1) < min(a, 1)) {
        counter    <- counter + 1
        dmat[j, ]  <- xs
        st         <- sn
      }
    }
    accept_rate[i] <- counter / nobs
    if (i %% 100 == 0) prg_bar()
  }
  if (warmup > 0) {
    theta_mat <- theta_mat[-c(1:warmup), , drop = FALSE]
    accept_rate <- accept_rate[-c(1:warmup)]
  }
  list(sample = theta_mat, accept_rate = accept_rate)
}



#' Summarise dpout object.
#'
#' @param gdpout object
#'
#' @return summary table
#' @export
summary.gdpout <- function(object) {
  posterior::summarise_draws(object$chain)
}

#' Summarise dpout object.
#'
#' @param dpout object
#'
#' @return trace plots
#' @export
plot.gdpout <- function(object) {
  bayesplot::mcmc_trace(object$chain)
}


