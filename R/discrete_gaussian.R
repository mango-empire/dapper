#' The Discrete Gaussian Distribution
#'
#' @param x vector of quantiles.
#' @param n number of random deviates.
#' @param mu location parameter.
#' @param sigma scale parameter.
#' @param log logical; if \code{TRUE}, probabilities are given as log(p).
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' P[X = x] = \dfrac{e^{-(x - \mu)^2/2\sigma^2}}{\sum_{y \in \mathbb{Z}} e^{-(x-\mu)^2/2\sigma^2}}.
#' }
#'
#' @references
#' Canonne, C. L., Kamath, G., & Steinke, T. (2020). The Discrete Gaussian for Differential Privacy.
#' \emph{arXiv}. <https://doi.org/10.48550/ARXIV.2004.00010>
#'
#' @return ddnorm gives the probability mass function and rdnorm
#' generates random deviates.
#'
#' @name DiscreteGuassian
#' @aliases DiscreteGuassian
#' @export
ddnorm <- function(x, mu = 0, sigma = 1, log = FALSE) {
    t1 <- exp(-(x - mu)^2 / (2 * sigma^2))
    t2 <- ddnorm_constant(sigma)
    if(log) {
      log (t1) - log (t2)
    } else {
      t1 / t2
    }
}

ddnorm_constant <- function(sigma) {
    fsum <- NULL
    psum <- NULL
    if(sigma^2 <= 1) {
        fsum <- 0
        t1 <- sum(sapply(1:1000, function(s) exp(-s^2 / (2 * sigma^2))))
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
        y <- dapper::rdlaplace(1, scale = t)
        c <- stats::rbinom(1, 1, exp(-(abs(y) - sigma^2/t)^2/(2*sigma^2)))
        if(c == 1) {
          smp[i] <- y
          break
        }
      }
    }
    smp + mu
}

