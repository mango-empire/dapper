#reproduce probit example from paper

set.seed(1)
theta <- -3.7
n <- 10^4
p <- pnorm(-3.7)
z <- rnorm(n, theta, 1)
y <- ifelse(z > 0, 1, 0)
r0 <- 1
b0 <- 0
niter <- 600

theta_vec <- numeric(niter)

gen_z <- function(y, t) {
  if(y == 1) {
    truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = t + b0, sd = sqrt(r0))
  } else {
    truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = t + b0, sd = sqrt(r0))
  }
}

tcur <- theta_vec[1]
for(i in 1:niter) {
  theta_vec[i] <- tcur
  z <- sapply(y, function(s) gen_z(s, tcur))
  tcur <- rnorm(1 , n^(-1) * sum(z - b0), sd = sqrt(n^(-1) * r0))
}

