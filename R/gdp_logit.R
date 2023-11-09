gdp_logit <- function(formula, data, pf, mu0, tau0) {
  tmp  <- parse_formula(formula, data)
  xmat <- tmp[[2]] #design matrix

  lik_factory <- function(xmat) {
    function(theta) {
      eta <- xmat %*% theta
      p   <- 1 / (1 + exp(-eta))
      matrix(rbinom(nrow(xmat), 1, p), ncol = 1)
    }
  }

  st_calc <- function(dmat) {
    c(dmat)
  }

  ll_priv_mech <- function(sdp, x) {
    t1 <- sum(sdp == x)
    t2 <- sum(sdp != x)
    t1 * log(pf + (1-pf) * 1/2) + t2 * log((1-pf) * 1/2)
  }

  logpost <- function(p, y, x) {
    eta <- x %*% p
    pp <- 1/(1 + exp(-eta))
    lp <- sum(dbinom(y, 1, pp, log= TRUE))
    lp <- lp + sum(dnorm(p, mu0, tau0, log = TRUE))
    lp
  }


  post_factory <- function(xmat) {
    function(dmat, theta) {
      t1 <- fmcmc::MCMC(logpost, initial = theta,
                        kernel = fmcmc::kernel_adapt(),
                        nsteps = 50,
                        burnin = 49,
                        progress = FALSE,
                        y = c(dmat),
                        x = xmat)

      t1[1,]
    }
  }

  new_privacy(post_smpl = post_factory(xmat),
              lik_smpl = lik_factory(xmat),
              ll_priv_mech = ll_priv_mech,
              st_calc = st_calc,
              add = FALSE,
              npar = ncol(xmat),
              varnames = tmp[[3]])

}
