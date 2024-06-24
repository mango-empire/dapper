#' Creates a data model.
#'
#' @param post_f a function that draws posterior samples given the confidential data.
#' @param latent_f a function that represents the latent data sampling model.
#' @param priv_f a function that represents the log likelihood of the privacy mechanism.
#' @param st_f a function that calculates the statistic to be released.
#' @param npar dimension of the parameter being estimated.
#' @param varnames an optional character vector of parameter names. Used to label summary outputs.
#'
#' @details
#' * `post_f()` is a function which makes draws from the posterior sampler. It has
#' the syntax `post_f(dmat, theta)`. Here `dmat` is an R matrix representing the confidential data.
#' Note `dmat` must be a matrix even if there is only one dimension. Thus, `dmat` cannot
#' be a vector for instance. This function can be constructed by wrapping MCMC samplers generated from other R packages
#' (e.g. \CRANpkg{rstan}, \CRANpkg{fmcmc}, \CRANpkg{adaptMCMC}).
#' If using this approach, it is recommended to avoid using packages
#' with a large initialization overhead such as \CRANpkg{mcmc} since the sampler is reinitialized
#' every loop iteration. In the case of \CRANpkg{mcmc},
#' the Metropolis-Hastings loop is implemented in C so there is a significant initialization cost
#' when calling from an R function. The `theta` argument is an R vector and its purpose is
#' to serve as the initialization point for `post_f`.
#'
#' * `priv_f()` is an R function that represents the log of the privacy mechanism density.
#' This function has the form `priv_f(sdp, sx)` where `sdp` and `sx` are both either
#' a R vector or matrix. The arguments must appear in the exact order with the same variables names as defined above.
#' Finally, the return value of `priv_f` must be a real number.
#'
#' * `st_f()` is an R function which calculates a summary statistic. It
#' has the syntax `st_f(i, xi, sdp)` where the three arguments must appear in the stated order.
#' The role of this function is to represent terms in the definition of record additivity.
#' Here the type class for `i` is an integer,
#' while `xi` is an R vector and `sdp` is an R vector or matrix.
#'
#' * `npar` is an integer value that represents the dimension of `theta`. It
#' should output a vector of the same length as `post_f`.
#' @md
#'
#' @return A S3 object of class \code{privacy}.
#' @export
#'
new_privacy <- function(post_f   = NULL,
                        latent_f = NULL,
                        priv_f   = NULL,
                        st_f     = NULL,
                        npar     = NULL,
                        varnames = NULL)
{
  checkmate::assert_function(post_f)
  checkmate::assert_function(latent_f)
  checkmate::assert_function(priv_f)
  checkmate::assert_function(st_f)
  checkmate::assert_count(npar)
  checkmate::assert(checkmate::check_character(varnames),
                    checkmate::qtest(varnames, "0"))


  plist <- list(post_f   = post_f,
                latent_f = latent_f,
                priv_f   = priv_f,
                st_f     = st_f,
                npar     = npar,
                varnames = varnames)
  structure(plist, class = "privacy")
}
