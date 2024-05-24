post_f <- function(dmat, theta) {
  n <- length(c(dmat))
  s <- sum(dmat)
  rbeta(1, shape1 = s + 1, shape2 = n - s + 1)
}

latent_f <- function(theta) {
  as.matrix(rbinom(4000, 1, theta), ncol = 1)
}

st_f <- function(i, xi, sdp) {
  x    <- rep(0, 4000)
  x[i] <- xi
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
x_conf <- rbinom(4000, 1, .5)
idp <- as.logical(rbinom(4000, 1, .5))

x_sdp <- x_conf
x_sdp[idp] <- rbinom(sum(idp), 1, .5)

with_progress({
out <- dapper_sample(dmod,
                     sdp = x_sdp,
                     init_par = .8,
                     niter = 500)
})

profvis({
dapper_chain(dmod,
              sdp = x_sdp,
              init_par = .8,
              niter = 10)

})

dmat <- latent_f(.25)
st_f(1, dmat)
reduce(lapply(1:nrow(dmat), function(i) st_f(i, dmat)), "+")
