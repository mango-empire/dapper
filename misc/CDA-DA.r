nsim <- 1000
eps <- 1
sig <- .01

f1 <- function(s, theta) {
  s2i <- eps^2 + 1/sig^2
  s2 <- 1/s2i
  mu <- (eps^2 * s + (1/sig^2) * theta)/(s2i)
  rnorm(1, mean = mu, sd = sqrt(s2))
}

f2 <- function(x) {
  rnorm(1, mean = x, sd = sig)
}

est_vec <- numeric(nsim)
est_vec[1] <- 1
s <- 10
for(i in 2:nsim) {
  x <- f1(s, est_vec[i-1])
  est_vec[i] <- f2(x)
}

plot(est_vec, type = "l")

acf(est_vec)

#----------------------------------------------------------

b0 <- 0
r0 <- sig^2 + eps^(-2)

accept <- function(to, tn, s) {
  k1 <- dnorm(s, mean = tn, sd = sqrt(sig^2 + eps^(-2)), log = TRUE)
  k2 <- dnorm(s, mean = to, sd = sqrt(sig^2 + eps^(-2)), log = TRUE)
  k3 <- dnorm(s, mean = tn, sd = sqrt(r0 + eps^(-2)), log = TRUE)
  k4 <- dnorm(s, mean = to, sd = sqrt(r0 + eps^(-2)), log = TRUE)
  log(runif(1)) < (k1 - k2 + k4 - k3)
}


est_vec <- numeric(nsim)
est_vec[1] <- 1
s <- 10
for(i in 2:nsim) {
  x <- f1(s, est_vec[i-1])
  tprop <- rnorm(1, mean = x - b0, sd = sqrt(r0))
  if(accept(est_vec[i-1], tprop, s)) {
    est_vec[i] <- tprop
  } else {
    est_vec[i] <- est_vec[i-1]
  }
}

plot(est_vec, type = "l")

acf(est_vec)










