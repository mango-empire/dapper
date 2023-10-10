#' Creates data model
#'
#' @param post_smpl A function that draws posterior samples.
#' @param lik_smpl The sampling model.
#' @param ll_priv_mech The log likelihood of the privacy mechanism.
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
                        prop_smpl = NULL,
                        ll_lik = NULL,
                        ll_prop = NULL,
                        ll_prior = NULL,
                        ll_priv_mech = NULL,
                        st_update = NULL,
                        st_calc = NULL,
                        npar = NULL)
{
  plist <- NULL
  if(is.null(prop_smpl)) {
    plist <- list(post_smpl = post_smpl,
                  lik_smpl = lik_smpl,
                  ll_priv_mech = ll_priv_mech,
                  st_update = st_update,
                  st_calc = st_calc,
                  npar = npar)
  } else {
    plist <- list(lik_smpl = lik_smpl,
                  prop_smpl = prop_smpl,
                  ll_lik = ll_lik,
                  ll_prop = ll_prop,
                  ll_prior = ll_prior,
                  ll_priv_mech = ll_priv_mech,
                  st_update = st_update,
                  st_calc = st_calc,
                  npar = npar)
  }
    structure(plist, class = "privacy")
}
