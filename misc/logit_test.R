library(MASS)
lm(bwt ~ smoke + race, data = birthwt)


data(birthwt)
posterior <- MCMClogit(low~age+as.factor(race)+smoke, data=birthwt)
plot(posterior)
summary(posterior)

posterior <- MCMClogit(low~smoke, data=birthwt)
plot(posterior)
summary(posterior)


tmp <- parse_formula(low~age+as.factor(race)+smoke, data=birthwt)

tmp <- gdp_logit(low~age+as.factor(race)+smoke, data=birthwt)

dmod <- gdp_logit(low~age+as.factor(race)+smoke, data=birthwt,pf = .70, mu0 = rep(0,5), tau0 = rep(20,5))


set.seed(1)
sdp <- birthwt$low
si  <- as.logical(rbinom(length(sdp),1,1- .70))
sdp[si] <- rbinom(sum(si),1,.5)

mean(sdp == birthwt$low)


library(progressr)
library(furrr)
plan(multisession, workers = 4)

with_progress({
tmp <- gdp_out <- gdp_sample(dmod,
                      sdp = sdp,
                      nobs = length(sdp),
                      niter = 20000,
                      warmup = 5000,
                      chains = 4,
                      init_par = rep(0,5))})




tmp <- gdp_out <- gdp_sample(dmod,
                             sdp = sdp,
                             nobs = length(sdp),
                             niter = 1000,
                             chains = 1,
                             init_par = rep(0,5))


 logpost <- function(p, y, x) {
  eta <- x %*% p
  pp <- 1/(1 + exp(-eta))
  lp <- sum(dbinom(y, 1, pp, log= TRUE))
  lp <- lp + sum(dnorm(p, mu0, tau0, log = TRUE))
  lp
}



tmp <- parse_formula(low~age+as.factor(race)+smoke, data=birthwt)
y    <- tmp[[1]] #response
x    <- tmp[[2]] #design matrix
dmat <- cbind(y,x)

mu0 <- rep(0, 5)
tau0 <- rep(20, 5)
tmp <- fmcmc::MCMC(logpost, initial = rep(0,5),
                    kernel = fmcmc::kernel_adapt(),
                    nsteps = 1000,
                    burnin = 10,
                    progress = FALSE,
                    y = c(dmat[,1]),
                    x = dmat[,-1])

plot(tmp)
summary(tmp)
logpost(rep(0,5), c(dmat[,1]), dmat[,-1])

