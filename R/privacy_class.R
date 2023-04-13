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
