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
#' the syntax `post_f(dmat, theta)`. Here `dmat` is a numeric matrix representing the confidential database
#' and `theta` is a numeric vector which serves as the initialization point for a one sample draw.
#' The easiest, bug-free way to construct `post_f()` is to use a conjugate prior. However,
#' this function can also be constructed by wrapping a MCMC sampler generated from other R packages
#' (e.g. \CRANpkg{rstan}, \CRANpkg{fmcmc}, \CRANpkg{adaptMCMC}).
#'
#' * `priv_f()` is a function that represents the log of the privacy mechanism density.
#' This function has the form `priv_f(sdp, sx)` where `sdp` and `sx` are both either
#' a numeric vector or matrix. The arguments must appear in the exact stated order with the same variables names as mentioned.
#' Finally, the return value of `priv_f()` must be a numeric vector of length one.
#'
#' * `st_f()` is a function which calculates a summary statistic. It
#' has the syntax `st_f(i, xi, sdp)` where the three arguments must appear in the stated order.
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

  assert_privacy <- checkmate::makeAssertionFunction(check_privacy)
  assert_privacy(list(post_f = post_f,
                     latent_f = latent_f,
                     st_f = st_f,
                     priv_f = priv_f))


  plist <- list(post_f   = post_f,
                latent_f = latent_f,
                priv_f   = priv_f,
                st_f     = st_f,
                npar     = npar,
                varnames = varnames)
  structure(plist, class = "privacy")
}

check_privacy <- function(x) {
    post_f <- x$post_f
    latent_f <- x$latent_f
    priv_f <- x$priv_f
    st_f <- x$st_f

    formal_args <- function(x) names(formals(x))

    if(!identical(formal_args(post_f), c("dmat", "theta"))) {
      return("Arguments of post_f() must be named (dmat) and (theta) in that order")
    }

    if(!identical(formal_args(st_f), c("xi", "sdp", "i"))) {
      return("Arguments of st_f() must be named (xi), (sdp) and (i) in that order")
    }

    if(!identical(formal_args(priv_f), c("sdp", "sx"))) {
      return("Arguments of priv_f() must be named (sdp) and (sx) in that order")
    }

    if(!identical(formal_args(latent_f), "theta")) {
      return("latent_f() must have an argument named (theta)")
    }

    TRUE
}
