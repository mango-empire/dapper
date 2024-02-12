#' Creates data model
#'
#' @param post_f An R function that draws posterior samples given the confidential data.
#' @param latent_f The latent data sampling model.
#' @param priv_f The log likelihood of the privacy mechanism modulo an additive constant.
#' @param st_f calculate statistic to be released.
#' @param add A logical argument. TRUE if st_f satisfies the record additivity assumption and FALSE otherwise.
#' @param npar Number of parameters in model.
#' @param varnames An optional character vector of parameter names. used to label summary outputs.
#'
#' @return A data model of class privacy. Is a S3 object.
#' @export
#'
new_privacy <- function(post_f   = NULL,
                        latent_f = NULL,
                        priv_f   = NULL,
                        st_f     = NULL,
                        add      = FALSE,
                        npar     = NULL,
                        varnames = NULL)
{
  checkmate::assert_function(post_f)
  checkmate::assert_function(latent_f)
  checkmate::assert_function(priv_f)
  checkmate::assert_function(st_f)
  checkmate::assert_logical(add)
  if(!is.null(st_f)) checkmate::assert_function(st_f)
  checkmate::assert_int(npar)

  plist <- list(post_f   = post_f,
                latent_f = latent_f,
                priv_f   = priv_f,
                st_f     = st_f,
                add      = add,
                npar     = npar,
                varnames = varnames)
  structure(plist, class = "privacy")
}
