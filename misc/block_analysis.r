#aggregate analysis
post_smpl <- function(dmat, theta) {
    B <- c(dmat)
    N <- theta[1]
    s1 <- B[1] + B[2]
    s2 <- B[3] + B[4]
    alpha <- rbeta(1, s1 + 1, s2 + 1)
    theta1 <- rbeta(1, B[1] + 1, s1 - B[1] + 1)
    theta2 <- rbeta(1, B[3] + 1, s2 - B[3] + 1)
    c(N, alpha, theta1, theta2)
}

lik_smpl <- function(theta) {
    N <- theta[1]
    alpha <- theta[2]
    theta1 <- theta[3]
    theta2 <- theta[4]
    s1 <- rbinom(1, N, alpha)
    s2 <- N - s1
    B <- numeric(4)

    B[1] <- rbinom(1, s1, theta1)
    B[2] <- s1 - B[1]
    B[3] <- rbinom(1, s2, theta2)
    B[4] <- s2 - B[3]

    B
}

gen_priv <- function(epsilon) {
    function(sdp, zt) {
        ep <- 1/epsilon
        -sum((sdp - zt)^2)/ep^2
    }
}


st_calc <- function(dmat) {
    dmat
}



set.seed(1)
alpha <- runif(1)
theta1 <- runif(1)
theta2 <- runif(1)

B <- numeric(4)
N <- 2000
s <- rbinom(1, N, alpha)
B[1] <- rbinom(1, s, theta1)
B[2] <- s - B[1]
B[3] <- rbinom(1, N-s, theta2)
B[4] <- N - s - B[3]


sigma <- 200
dp_noise <-  rnorm(4, mean = 0, sd = sigma)
dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = gen_priv(1/sigma),
                    st_update = NULL,
                    st_calc = st_calc,
                    npar = 4)

tmp <- mcmc_privacy(dmod,
                    sdp = B + dp_noise,
                    nobs = 1,
                    init_par = c(sum(B), .5, .5, .5),
                    niter = 100000,
                    chains = 1,
                    varnames = c('N', 'alpha', 'theta1', 'theta2'))

summary(tmp)


#individual data


post_smpl <- function(dmat, theta) {
    N <- theta[1]
    B <- apply(dmat, 2, sum)
    B <- c(B, N - sum(B))
    s1 <- B[1] + B[2]
    s2 <- B[3] + B[4]
    alpha <- rbeta(1, s1 + 1, s2 + 1)
    theta1 <- rbeta(1, B[1] + 1, s1 - B[1] + 1)
    theta2 <- rbeta(1, B[3] + 1, s2 - B[3] + 1)
    c(N, alpha, theta1, theta2)
}

lik_smpl <- function(theta) {
    N <- theta[1]
    alpha <- theta[2]
    theta1 <- theta[3]
    theta2 <- theta[4]
    s1 <- rbinom(1, 1, alpha)

    B <- numeric(3)
    t1 <- rbinom(1, 1, theta1)
    t2 <- rbinom(1, 1, theta2)
    B[1] <- as.integer(s1 & t1)
    B[2] <- as.integer(s1 & !t1)
    B[3] <- as.integer(!s1 & t2)

    B
}



gen_priv <- function(epsilon) {
    function(sdp, zt) {
        ep <- 1/epsilon
        -sum((sdp - zt)^2)/ep^2
    }
}


st_calc <- function(dmat) {
    N <- nrow(dmat)
    B <- apply(dmat,2, sum)
    c(B, N - sum(B))
}



set.seed(1)
alpha <- runif(1)
theta1 <- runif(1)
theta2 <- runif(1)

B <- numeric(4)
N <- 200
s <- rbinom(1, N, alpha)
B[1] <- rbinom(1, s, theta1)
B[2] <- s - B[1]
B[3] <- rbinom(1, N-s, theta2)
B[4] <- N - s - B[3]


sigma <- 40
dp_noise <-  rnorm(4, mean = 0, sd = sigma)

dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = gen_priv(1/sigma),
                    st_update = NULL,
                    st_calc = st_calc,
                    npar = 4)

tmp <- mcmc_privacy(dmod,
                    sdp = B + dp_noise,
                    nobs = N,
                    init_par = c(sum(B), .5, .5, .5),
                    niter = 5000,
                    chains = 1)

summary(tmp)

