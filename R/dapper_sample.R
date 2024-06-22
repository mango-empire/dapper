#' Private Posterior Sampler
#'
#' @description
#' Generates samples from the private posterior using a data augmentation framework.
#'
#' @param data_model a data model represented by a privacy class object.
#' @param sdp the observed privatized data. Must be a vector or matrix.
#' @param init_par initial starting point of the chain.
#' @param niter number of draws.
#' @param warmup number of iterations to discard as warmup. Default is half of niter.
#' @param chains number of MCMC chains to run. Can be done in parallel or sequentially.
#'
#' @return A dpout object which contains:
#' \item{accept_prob}{acceptance probabilities.}
#' \item{chain}{samples from the private posterior.}
#' @export
#'
#' @references
#' Ju, N., Awan, J. A., Gong, R., & Rao, V. A. (2022). Data Augmentation MCMC
#' for Bayesian Inference from Privatized Data. \emph{arXiv}. <https://doi.org/10.48550/ARXIV.2206.00710>
#'
#' @examples
#' #simulate confidential data
#' #privacy mechanism adds gaussian noise to each observation.
#' set.seed(1)
#' n <- 100
#' eps <- 3
#' y <- rnorm(n, mean = -2, sd = 1)
#' sdp <- mean(y) + rnorm(1, 0, 1/eps)
#'
#' post_f <- function(dmat, theta) {
#'     x <- c(dmat)
#'     xbar <- mean(x)
#'     n <- length(x)
#'     pr_m <- 0
#'     pr_s2 <- 4
#'     ps_s2 <- 1/(1/pr_s2 + n)
#'     ps_m <- ps_s2 * ((1/pr_s2)*pr_m + n * xbar)
#'     rnorm(1, mean = ps_m, sd = sqrt(ps_s2))
#' }
#' latent_f <- function(theta) {
#'     matrix(rnorm(100, mean = theta, sd = 1), ncol = 1)
#' }
#' st_f <- function(xi, sdp, i) {
#'     mean(xi)
#' }
#' priv_f <- function(sdp, tx) {
#'   sum(dnorm(sdp - tx, 0, 1/eps, TRUE))
#' }
#' dmod <- new_privacy(post_f = post_f,
#'   latent_f = latent_f,
#'   priv_f = priv_f,
#'   st_f = st_f,
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
  checkmate::assert_class(data_model, "privacy")
  checkmate::qassert(chains, "X?(0,)")
  if(length(init_par) != data_model$npar) stop("Dimension of initial parameter does not match privacy model")


  p <- progressr::progressor(niter * chains)
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
  st          <- Reduce("+", lapply(1:nrow(dmat), function(i) st_f(dmat[i,], sdp, i)))
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

      sn <- st - st_f(xo, sdp, j) + st_f(xs, sdp, j)
      a <- priv_f(sdp, sn) - priv_f(sdp, st)
      if (log(stats::runif(1)) < a) {
        counter    <- counter + 1
        dmat[j, ]  <- xs
        st         <- sn
      }
    }
    accept_rate[i] <- counter / nobs
    if (!is.null(prg_bar)) prg_bar(message = sprintf("Iteration %g", i))
  }

  if (warmup > 0) {
    theta_mat <- theta_mat[-c(1:warmup), , drop = FALSE]
    accept_rate <- accept_rate[-c(1:warmup)]
  }

  list(sample = theta_mat, accept_rate = accept_rate)
}



#' Summarise dpout object.
#'
#' @param object dp_out object
#' @param ... optional arguments to `summarise_draws()`.
#'
#' @return a summary table of MCMC statistics.
#' @export
summary.dpout <- function(object, ...) {
  posterior::summarise_draws(object$chain, ...)
}

#' Plot dpout object.
#'
#' @param x dp_out object.
#' @param ... optional arguments to `mcmc_trace()`.
#' @return trace plots.
#' @export
plot.dpout <- function(x, ...) {
  bayesplot::mcmc_trace(x$chain, ...)
}


