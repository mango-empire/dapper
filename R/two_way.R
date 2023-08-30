#' DP aware inference for two way table.
#'
#' @param cell_values
#' @param n_total
#' @param init_par
#'
#' @return p1,p2
#' @export
#'
#' @examples
#' set.seed(1)
#' rho <- 3.325 * 3.01 * 10^(-2)
#' tval <- stats::rmultinom(1, 300, c(.20,.15,.35,.30))
#' sdp <- sapply(tval, function(s) rdpnorm(1, 0, sqrt(1/rho))) + tval
#' sout <- two_way(sdp, sum(tval), rho, 8000)
#' plot(sout)
two_way <- function(cell_values, n_total, rho, niter = 100, chains =  1) {
  dmod <- new_privacy(post_smpl = post_smpl_two_way,
                      lik_smpl = lik_smpl_two_way,
                      ll_priv_mech = gen_priv_two_way(sqrt(1/rho)),
                      st_update = NULL,
                      st_calc = st_calc_two_way,
                      npar = 3)

  sim_out <- mcmc_privacy(dmod,
               sdp = cell_values,
               nobs = 1,
               init_par = c(.5, .5, n_total),
               niter = niter,
               chains = chains)
  sim_out
}

# table layout:
# a b
# c d
# dmat is a row vector: [a,b,c,d]
# theta format: [p1, p2, total_count]
post_smpl_two_way <- function(dmat, theta) {
  x <- c(dmat)
  n <- theta[3]
  alpha <- rep(1/2,4)

  k1 <- x[1] + x[3]
  a1 <- k1 + alpha[1]
  b1 <- n - k1 + alpha[2]
  p1 <- stats::rbeta(1, a1, b1)

  k2 <- x[1] + x[2]
  a2 <- k2 + alpha[3]
  b2 <- n - k2 + alpha[4]
  p2 <- stats::rbeta(1, a2, b2)

  c(p1, p2, n)
}

lik_smpl_two_way <- function(theta) {
  p11 <- theta[1] * theta[2]
  p12 <- (1 - theta[1]) * theta[2]
  p21 <- theta[1] * (1 - theta[2])
  p22 <- (1-theta[1]) * (1 - theta[2])
  n <- theta[3]
  c(stats::rmultinom(1, n, c(p11, p12, p21, p22)))
}

gen_priv_two_way <- function(sd) {
  function(sdp, zt) {
    sum(-(sdp - zt)^2/(2 * sd^2))
  }
}

st_calc_two_way <- function(dmat) {
  c(dmat)
}



