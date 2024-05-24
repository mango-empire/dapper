#' Generate samples from the private posterior.
#'
#' @param data_model A data model represented by a privacy class object.
#' @param sdp The observed privatized data.
#' @param init_par Initial starting point of the chain.
#' @param niter Number of draws.
#' @param warmup Number of iterations to discard as warmup. Default is half of niter.
#'
#' @return A dpout object which contains:
#' \item{accept_prob}{Acceptance probabilities}
#' \item{chain}{Samples from the private posterior}
#' @export
#'
#' @examples
#' #simulate confidential data
#' #Privacy mechanism adds gaussian noise to each observation.
#' set.seed(1)
#' n <- 100
#' eps <- 3
#' y <- rnorm(n, mean = -2, sd = 1)
#' sdp <- mean(y) + rnorm(1, 0, 1/eps)
#'
#' post_smpl <- function(dmat, theta) {
#'     x <- c(dmat)
#'     xbar <- mean(x)
#'     n <- length(x)
#'     pr_m <- 0
#'     pr_s2 <- 4
#'     ps_s2 <- 1/(1/pr_s2 + n)
#'     ps_m <- ps_s2 * ((1/pr_s2)*pr_m + n * xbar)
#'     rnorm(1, mean = ps_m, sd = sqrt(ps_s2))
#' }
#' lik_smpl <- function(theta) {
#'     matrix(rnorm(100, mean = theta, sd = 1), ncol = 1)
#' }
#' st_calc <- function(dmat) {
#'     mean(dmat)
#' }
#' priv_mech <- function(sdp, zt) {
#'   sum(dnorm(sdp - zt, 0, 1/eps, TRUE))
#' }
#' dmod <- new_privacy(post_f = post_smpl,
#'   latent_f = lik_smpl,
#'   priv_f = priv_mech,
#'   st_f = st_calc,
#'   add = FALSE,
#'   npar = 1)
#'
#' out <- dapper_sample(dmod,
#'                     sdp = sdp,
#'                     init_par = -2,
#'                     niter = 500)
#' summary(out)
dapper_sample <- function(data_model,
                       sdp,
                       init_par,
                       niter = 2000,
                       warmup = floor(niter / 2),
                       chains = 1) {
  #check inputs
  checkmate::qassert(chains, "X?(0,)")
  if(length(init_par) != data_model$npar) stop("Dimension of initial parameter does not match privacy model")

  pb_size <- floor(niter / 100)
  p <- progressr::progressor(pb_size * chains)
  fout <- furrr::future_map(rep(niter, chains), dapper_chain,
                            data_model = data_model,
                            sdp = sdp,
                            init_par = init_par,
                            warmup = warmup,
                            prg_bar = p,
                            .options = furrr::furrr_options(seed = TRUE))

  theta_clist <- lapply(1:chains, function(s) fout[[s]]$sample)
  accept_mat <- do.call(cbind, lapply(1:chains, function(s) fout[[s]]$accept_rate))
  new_dpout(theta_clist, accept_mat, data_model$varnames)
}


#' single chain sample
#'
#' @export
dapper_chain <- function(data_model,
                      sdp,
                      init_par,
                      niter = 2000,
                      warmup = floor(niter / 2),
                      prg_bar = NULL) {
  #check inputs
  checkmate::qassert(niter, "X?(0,)")
  checkmate::assert_class(data_model, "privacy")
  post_f    <- data_model$post_f
  latent_f  <- data_model$latent_f
  priv_f    <- data_model$priv_f
  st_f      <- data_model$st_f
  npar      <- data_model$npar

  accept_rate <- numeric(niter)
  theta_clist <- list()
  dmat        <- latent_f(init_par)
  theta_mat   <- matrix(0, nrow = niter, ncol = npar)
  theta       <- init_par
  st          <- apply(sapply(1:nrow(dmat), function(i) st_f(i, dmat[i,], sdp)), 1, sum)
  nobs        <- nrow(dmat)

  for (i in 1:niter) {
    counter        <- 0
    theta_mat[i, ] <- theta
    theta          <- post_f(dmat, theta)
    smat           <- latent_f(theta)

    for (j in 1:nobs) {
      xs <- smat[j, ]
      xo <- dmat[j, ]
      sn <- NULL

      sn <- st - st_f(j, xo, sdp) + st_f(j, xs, sdp)

      a <- exp(priv_f(sdp, sn) - priv_f(sdp, st))
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
#' @param dpout object
#'
#' @return summary table
#' @export
summary.dpout <- function(object) {
  posterior::summarise_draws(object$chain)
}

#' Summarise dpout object.
#'
#' @param dpout object
#'
#' @return trace plots
#' @export
plot.dpout <- function(object) {
  bayesplot::mcmc_trace(object$chain)
}


