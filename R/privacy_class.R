#' Creates data model
#'
#' @param post_f a function that draws posterior samples given the confidential data.
#' @param latent_f a function that represents the latent data sampling model.
#' @param priv_f a function that represents the log likelihood of the privacy mechanism.
#' @param st_f a function that calculated the statistic to be released.
#' @param npar number of parameters in model.
#' @param varnames an optional character vector of parameter names. Used to label summary outputs.
#'
#' @return A data model of class privacy. Is a S3 object.
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
  if(!is.null(st_f)) checkmate::assert_function(st_f)
  checkmate::assert_int(npar)

  plist <- list(post_f   = post_f,
                latent_f = latent_f,
                priv_f   = priv_f,
                st_f     = st_f,
                npar     = npar,
                varnames = varnames)
  structure(plist, class = "privacy")
}
