#' Discrete Laplace Distribution
#'
#' @param x a vector of quantiles.
#' @param scale the scale parameter.
#' @param log logical; if TRUE, probabilities are given as log(p).
#'
#' @return dnorm gives the probability mass function and rdnorm
#' generates random deviates.
#' @export
ddlaplace <- function(x, scale = 1, log = FALSE) {
    s  <- scale
    t1 <- log(exp(1/s) - 1) - log(exp(1/s) + 1) - abs(x)/s
    if(log) {
        t1
    } else {
        exp(t1)
    }
}

rdlaplace <- function(n, scale = 1) {
    t <- scale
    smp <- numeric(n)
    for(i in 1:n) {
      while(TRUE) {
          while(TRUE) {
            u <- sample(0:(t-1), 1)
            d <- rbinom(1, 1, exp(-u/t))
            if(d) break
          }
          v <- 0
          while(TRUE) {
            a <- rbinom(1, 1, exp(-1))
            if(!a) {
              break
            } else {
              v <- v + 1
            }
          }
          y <- u + t * v
          b <- rbinom(1, 1, 1/2)
          if(b == 0 | y != 0) {
            smp[i] <- (1 - 2 * b) * y
            break
          }
      }
    }
    smp
}


