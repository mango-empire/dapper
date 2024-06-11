#' The Discrete Gaussian Distribution
#'
#' @param x vector of quantiles.
#' @param mu location parameter
#' @param sigma scale parameter
#' @param log logical; if TRUE, probabilities are given as log(p).
#'
#' @return dnorm gives the probability mass function and rdnorm
#' generates random deviates.
#' @export
ddnorm <- function(x, mu = 0, sigma = 1, log = FALSE) {
    si <- 10 * min(10, sigma)
    t1 <- exp(-(x - mu)^2 / (2 * sigma^2))
    t2 <- sapply(-si:si, function(s) exp(-(s - mu)^2 / (2 * sigma^2)))
    if(log) {
      log (t1 / sum (t2))
    } else {
      t1 / sum(t2)
    }
}

rdnorm <- function(n, mu = 0, sigma = 1) {
    t <- floor(sigma) + 1
    smp <- numeric(n)
    for(i in 1:n) {
      while(TRUE) {
        y <- rdlaplace(1, scale = t)
        c <- rbinom(1, 1, exp(-(abs(y) - sigma^2/t)^2/(2*sigma^2)))
        if(c == 1) {
          smp[i] <- y
          break
        }
      }
    }
    smp
}

