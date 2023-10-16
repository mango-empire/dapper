#' Creates data model
#'
#' @param post_smpl A function that draws posterior samples.
#' @param lik_smpl The sampling model.
#' @param ll_priv_mech The log likelihood of the privacy mechanism modulo an additive constant.
#' @param st_update Update function for statistic.
#' @param st_calc calculate st statistic.
#' @param npar Number of parameters.
#'
#' @return A data model of class privacy. Is a S3 object.
#' @export
#'
#' @examples
new_privacy <- function(post_smpl = NULL,
                        lik_smpl = NULL,
                        ll_priv_mech = NULL,
                        st_calc = NULL,
                        add = FALSE,
                        npar = NULL)
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
                  npar = npar)
  structure(plist, class = "privacy")
}
