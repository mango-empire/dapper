#' Creates data model
#'
#' @param post_smpl A function that draws posterior samples given the confidential data.
#' @param lik_smpl The sampling model.
#' @param ll_priv_mech The log likelihood of the privacy mechanism modulo an additive constant.
#' @param st_calc calculate statistic to be released.
#' @param npar Number of parameters in model.
#'
#' @return A data model of class privacy. Is a S3 object.
#' @export
#'
new_privacy <- function(post_smpl = NULL,
                        lik_smpl = NULL,
                        ll_priv_mech = NULL,
                        st_calc = NULL,
                        add = FALSE,
                        npar = NULL,
                        varnames = NULL)
{
  checkmate::assert_function(post_smpl)
  checkmate::assert_function(lik_smpl)
  checkmate::assert_function(ll_priv_mech)
  checkmate::assert_function(st_calc)
  checkmate::assert_logical(add)
  if(!is.null(st_calc)) checkmate::assert_function(st_calc)
  checkmate::assert_int(npar)

  plist <- list(post_smpl = post_smpl,
                  lik_smpl = lik_smpl,
                  ll_priv_mech = ll_priv_mech,
                  st_calc = st_calc,
                  add = add,
                  npar = npar,
                  varnames = varnames)
  structure(plist, class = "privacy")
}
