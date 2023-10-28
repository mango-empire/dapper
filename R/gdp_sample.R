
#' Sample from posterior given private data
#'
#' @param data_model A data model created using the privacy class.
#' @param sdp The observed privatized data.
#' @param nobs The number of observations.
#' @param init_par Initial starting point of the chain.
#' @param niter Number of draws.
#' @param warmup Number of iterations to discard as warmup. Default is half of niter.
#' @param varnames Optional character vector specifying parameter names in the output.
#'
#' @return A matrix of posterior samples.
#' @export
#'
#' @examples
gdp_sample <- function(data_model, sdp, nobs, init_par,
                       niter = 2000, warmup = floor(niter / 2), varnames = NULL) {
  #check inputs
  checkmate::qassert(nobs, "X?(0,)")
  checkmate::qassert(niter, "X?(0,)")
  checkmate::qassert(varnames, c("0", "s+"))
  checkmate::assert_class(data_model, "privacy")
  if(length(init_par) != data_model$npar) stop("Dimension of initial parameter does not match privacy model")


  post_smpl <- data_model$post_smpl
  lik_smpl <- data_model$lik_smpl
  ll_priv_mech <- data_model$ll_priv_mech
  st_calc <- data_model$st_calc
  npar <- data_model$npar

  accept_rate <- numeric(niter)
  theta_clist <- list()
  ret_val <- NULL
  d_mat <- lapply(1:nobs, function(s)
  lik_smpl(init_par))
  d_mat <- do.call(rbind, d_mat)
  theta_mat <- matrix(0, nrow = niter, ncol = npar)
  theta <- init_par
  pb_size <- floor(niter / 100)
  pb <- progressr::progressor(pb_size)
  st <- st_calc(d_mat)
  for (i in 1:niter) {
    counter <- 0
    theta_mat[i, ] <- theta
    theta <- post_smpl(d_mat, theta)
    for (j in 1:nobs) {
      xs <- lik_smpl(theta)
      xo <- d_mat[j, ]
      sn <- NULL
      if (!data_model$add) {
        n_mat <- d_mat
        n_mat[j, ] <- xs
        sn <- st_calc(n_mat)
      } else {
        #sn <- st_update(st, xs, xo)
        sn <- st - st_calc(t(xo)) + st_calc(t(xs))
      }
      a <- exp(ll_priv_mech(sdp, sn) - ll_priv_mech(sdp, st))
      if (stats::runif(1) < min(a, 1)) {
        counter <- counter + 1
        d_mat[j, ] <- xs
        st <- sn
      }
    }
    accept_rate <- counter / nobs
    if (i %% 100 == 0)
      pb()
  }
  if (warmup > 0) {
    theta_clist[[1]] <- theta_mat[-c(1:warmup), , drop = FALSE]
    accept_rate <- accept_rate[-c(1:warmup)]
  }
  else {
    theta_clist[[1]] <- theta_mat
  }

  new_gdpout(theta_clist, accept_rate, varnames)
}

#' Summarise dpout object.
#'
#' @param dpout object
#'
#' @return summary table
#' @export
#'
#' @examples
summary.gdpout <- function(object) {
  print(paste0("Average Acceptance Probability: ", mean(object$accept_prob)))
  posterior::summarise_draws(object$chain)
}
