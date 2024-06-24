new_normal_privacy <- function() {
    post_f <- function(dmat, theta) {
      x <- c(dmat)
      xbar <- mean(x)
      n <- length(x)
      pr_m <- 0
      pr_s2 <- 4
      ps_s2 <- 1/(1/pr_s2 + n)
      ps_m <- ps_s2 * ((1/pr_s2)*pr_m + n * xbar)
      rnorm(1, mean = ps_m, sd = sqrt(ps_s2))
    }
    latent_f <- function(theta) {
      matrix(rnorm(100, mean = theta, sd = 1), ncol = 1)
    }
    st_f <- function(xi, sdp, i) {
      mean(xi)
    }
    priv_f <- function(sdp, tx) {
      sum(dnorm(sdp - tx, 0, 1/eps, TRUE))
    }
    list(post_f = post_f,
         latent_f = latent_f,
         st_f = st_f,
         priv_f = priv_f)
}
