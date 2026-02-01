#' `privacy` Object Constructor.
#'
#' @description
#' Creates a `privacy` object to be used as input into dapper_sample().
#'
#'
#' @param posterior_f a function that draws posterior samples given the confidential data.
#' @param latent_f a function that represents the latent data sampling model.
#' @param mechanism_f a function that represents the log likelihood of the privacy mechanism.
#' @param statistic_f a function that calculates the statistic to be released.
#' @param npar dimension of the parameter being estimated.
#' @param varnames an optional character vector of parameter names. Used to label summary outputs.
#'
#' @details
#' * posterior_f() is a function which makes draws from the posterior sampler. It has
#' the syntax posterior_f(dmat, theta). Here `dmat` is a numeric matrix representing the confidential database
#' and `theta` is a numeric vector which serves as the initialization point for a one sample draw.
#' The easiest, bug-free way to construct posterior_f() is to use a conjugate prior. However,
#' this function can also be constructed by wrapping a MCMC sampler generated from other R packages
#' (e.g. \CRANpkg{rstan}, \CRANpkg{fmcmc}, \CRANpkg{adaptMCMC}).
#'
#' * mechanism_f() is a function that represents the log of the privacy mechanism density.
#' This function has the form mechanism_f(sdp, sx) where `sdp` and `sx` are both either
#' a numeric vector or matrix. The arguments must appear in the exact stated order with the same variables names as mentioned.
#' Finally, the return value of mechanism_f() must be a numeric vector of length one.
#'
#' * statistic_f() is a function which calculates a summary statistic. It
#' has the syntax statistic_f(xi, sdp, i) where the three arguments must appear in the stated order.
#' The role of this function is to represent terms in the definition of record additivity.
#' Here `i` is an integer,
#' while `xi` is an numeric vector and `sdp` is a numeric vector or matrix.
#'
#' * `npar` is an integer equal to the dimension of `theta`.
#' @md
#'
#' @return A S3 object of class \code{privacy}.
#' @export
#'
new_privacy <- function(posterior_f = NULL,
                        latent_f    = NULL,
                        mechanism_f = NULL,
                        statistic_f = NULL,
                        npar        = NULL,
                        varnames    = NULL)
{
  checkmate::assert_function(posterior_f)
  checkmate::assert_function(latent_f)
  checkmate::assert_function(mechanism_f)
  checkmate::assert_function(statistic_f)
  checkmate::assert_count(npar)
  checkmate::assert(checkmate::check_character(varnames),
                    checkmate::qtest(varnames, "0"))

  assert_privacy <- checkmate::makeAssertionFunction(check_privacy)
  assert_privacy(list(posterior_f   = posterior_f,
                      latent_f      = latent_f,
                      statistic_f   = statistic_f,
                      mechanism_f   = mechanism_f))


  plist <- list(posterior_f = posterior_f,
                latent_f    = latent_f,
                mechanism_f = mechanism_f,
                statistic_f = statistic_f,
                npar        = npar,
                varnames    = varnames)
  structure(plist, class = "privacy")
}

check_privacy <- function(x) {
  
    post_f   <- x$posterior_f
    latent_f <- x$latent_f
    mech_f   <- x$mechanism_f
    stat_f   <- x$statistic_f

    formal_args <- function(x) names(formals(x))

    if(!identical(formal_args(post_f), c("dmat", "theta"))) {
      return("Arguments of posterior_f() must be named (dmat) and (theta) in that order")
    }

    if(!identical(formal_args(stat_f), c("xi", "sdp", "i"))) {
      return("Arguments of statistic_f() must be named (xi), (sdp) and (i) in that order")
    }

    if(!identical(formal_args(mech_f), c("sdp", "sx"))) {
      return("Arguments of mechanism_f() must be named (sdp) and (sx) in that order")
    }

    if(!identical(formal_args(latent_f), "theta")) {
      return("latent_f() must have an argument named (theta)")
    }

    TRUE
}
