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
two_way <- function(cell_values, n_total, epsilon, niter = 100) {
  dmod <- new_privacy(post_smpl = post_smpl_two_way,
                      lik_smpl = lik_smpl_two_way,
                      ll_priv_mech = gen_priv_two_way(epsilon),
                      st_update = NULL,
                      st_calc = st_calc_two_way,
                      npar = 3)

  sim_out <- mcmc_privacy(dmod,
               sdp = cell_values,
               nobs = 1,
               init_par = c(.5, .5, n_total),
               niter = niter)
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
  p1 <- rbeta(1, a1, b1)

  k2 <- x[1] + x[2]
  a2 <- k2 + alpha[3]
  b2 <- n - k2 + alpha[4]
  p2 <- rbeta(1, a2, b2)

  c(p1, p2, n)
}

lik_smpl_two_way <- function(theta) {
  p11 <- theta[1] * theta[2]
  p12 <- (1 - theta[1]) * theta[2]
  p21 <- theta[1] * (1 - theta[2])
  p22 <- (1-theta[1]) * (1 - theta[2])
  n <- theta[3]
  c(rmultinom(1, n, c(p11, p12, p21, p22)))
}

gen_priv_two_way <- function(epsilon) {
  function(sdp, zt) {
    ep <- 1/epsilon
    sum(-abs(sdp - zt)/ep)
  }
}

st_calc_two_way <- function(dmat) {
  c(dmat)
}



