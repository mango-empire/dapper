#' Creates data model
#'
#' @param post_smpl A function that draws posterior samples.
#' @param lik_smpl The sampling model.
#' @param ll_priv_mech The log likelihood of the privacy mechanism.
#' @param st_update Update function for statistic.
#' @param st_init Starting value of st.
#' @param npar Number of parameters.
#'
#' @return A data model of class privacy. Is a S3 object.
#' @export
#'
#' @examples ...
new_privacy <- function(post_smpl,
                        lik_smpl,
                        ll_priv_mech,
                        st_update,
                        st_init,
                        npar)
{
    plist <- list(post_smpl = post_smpl,
                  lik_smpl = lik_smpl,
                  ll_priv_mech = ll_priv_mech,
                  st_update = st_update,
                  st_init = st_init,
                  npar = npar)
    structure(plist, class = "privacy")
}
