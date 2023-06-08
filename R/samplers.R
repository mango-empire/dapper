rdplap <- function(n, t) {
  p <- exp(-1/t)
  DiscreteLaplace::rdlaplace(n, p, p)
}

rdpnorm_one <- function(s) {
  t <- floor(s) + 1
  while(TRUE) {
    y <- rdplap(1, t)
    p <- exp(-(abs(y) - s^2/t)^2/(2*s^2))
    b <- rbinom(1,1,p)
    if(b) return(y)
  }
}

rdpnorm <- function(n, mean = 0, sd = 1) {
  sapply(1:n, function(x) rdpnorm_one(sd)) + mean
}

ddpnorm <- function(x, mean = 0, sd = 1) {
  lb <- -floor(sd * 40)
  ub <- floor(sd * 40)
  K <- sum(sapply(lb:ub, function(s) exp(-(s - mean)^2/(2*sd^2))))
  exp(-(x - mean)^2/(2*sd^2))/K
}
