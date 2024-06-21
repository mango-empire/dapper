#' Discrete Laplace Distribution
#'
#' @param x a vector of quantiles.
#' @param n number of random deviates.
#' @param scale the scale parameter.
#' @param log logical; if \code{TRUE}, probabilities are given as log(p).
#'
#' @details
#' Probability mass function
#' \deqn{
#' P[X=x] = \dfrac{e^{1/t} - 1}{e^{1/t} + 1} e^{-|x|/t}.
#' }
#'
#' @examples
#' # mass function
#' ddlaplace(0)
#'
#' # mass function is vectorized
#' ddlaplace(0:10, scale = 5)
#'
#' # generate random samples
#' rdlaplace(10)
#'
#' @references
#' Canonne, C. L., Kamath, G., & Steinke, T. (2020). The Discrete Gaussian for Differential Privacy.
#' \emph{arXiv}. <https://doi.org/10.48550/ARXIV.2004.00010>
#'
#' @return ddlaplace gives the probability mass function and rdlaplace
#' generates random deviates.
#'
#' @name DiscreteLaplace
#'
#' @export

ddlaplace <- function(x, scale = 1, log = FALSE) {
    #check inputs
    checkmate::assertNumeric(x)
    checkmate::qassert(scale, "n1[0,)")
    checkmate::qassert(log, "b1")

    s  <- scale
    t1 <- log(exp(1/s) - 1) - log(exp(1/s) + 1) - abs(x)/s
    if(log) {
        t1
    } else {
        exp(t1)
    }
}

#' @rdname DiscreteLaplace
#' @export

rdlaplace <- function(n, scale = 1) {
    #check inputs
    checkmate::assertCount(n)
    checkmate::qassert(scale, "n1[0,)")

    t <- scale
    smp <- numeric(n)
    for(i in 1:n) {
      while(TRUE) {
          while(TRUE) {
            u <- sample(0:(t-1), 1)
            d <- stats::rbinom(1, 1, exp(-u/t))
            if(d) break
          }
          v <- 0
          while(TRUE) {
            a <- stats::rbinom(1, 1, exp(-1))
            if(!a) {
              break
            } else {
              v <- v + 1
            }
          }
          y <- u + t * v
          b <- stats::rbinom(1, 1, 1/2)
          if(b == 0 | y != 0) {
            smp[i] <- (1 - 2 * b) * y
            break
          }
      }
    }
    smp
}


