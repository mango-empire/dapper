post_f <- function(dmat, theta) {
  n <- length(c(dmat))
  s <- sum(dmat)
  rbeta(1, shape1 = s + 1, shape2 = n - 2 + 1)
}

latent_f <- function(theta) {
  t(rbinom(100, 1, theta))
}

st_f <- function(i, dmat) {
  t1 <- c(dmat)
  x    <- rep(0, 100)
  x[i] <- t1[i]
  x
}

priv_f <- function(sdp, lx) {
  t1 <- sum(sdp == lx)
  t2 <- sum(sdp != lx)
  t1 * log(3/4) + t2 * log(1/4)
}


dmod <- new_privacy(post_f = post_f,
                    latent_f = latent_f,
                    priv_f = priv_f,
                    st_f = st_f,
                    npar = 1)


set.seed(1)
x_conf <- rbinom(100, 1, .5)
idp <- as.logical(rbinom(100, 1, .5))

x_sdp <- x_conf
x_sdp[idp] <- rbinom(sum(idp), 1, .5)

out <- dapper_sample(dmod,
                     sdp = x_sdp,
                     init_par = .8,
                     niter = 5000)
