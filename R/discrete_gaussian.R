#' The Discrete Gaussian Distribution
#'
#' @param x vector of quantiles.
#' @param n number of random deviates.
#' @param mu location parameter
#' @param sigma scale parameter
#' @param log logical; if TRUE, probabilities are given as log(p).
#'
#' @return dnorm gives the probability mass function and rdnorm
#' generates random deviates.
#'
#' @name DiscreteGuassian
#' @aliases DiscreteGuassian
#' @export
#'
ddnorm <- function(x, mu = 0, sigma = 1, log = FALSE) {
    t1 <- exp(-(x - mu)^2 / (2 * sigma^2))
    t2 <- ddnorm_constant(mu, sigma)
    if(log) {
      log (t1) - log (t2)
    } else {
      t1 / t2
    }
}

ddnorm_constant <- function(mu, sigma) {
    fsum <- NULL
    psum <- NULL
    if(sigma^2 <= 1) {
        fsum <- 0
        t1 <- sum(sapply(1:1000, function(s) exp(-(s - mu)^2 / (2 * sigma^2))))
        fsum <- 2 * t1 + 1
    } else if(sigma^2 * 100 >= 1) {
        psum <- 0
        t2 <- sum(sapply(1:1000, function(s) exp(-pi^2 * sigma^2 * 2 * s^2)))
        psum <- sqrt(2 * pi * sigma^2) * (1 + 2 * t2)
    }

    if(is.null(fsum)) {
        psum
    } else if(is.null(psum)) {
        fsum
    } else {
        (psum + fsum) / 2
    }
}

#' @rdname DiscreteGuassian
#' @export

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

