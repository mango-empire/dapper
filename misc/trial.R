post_f <- function(dmat, theta) {
  n <- length(c(dmat))
  s <- sum(dmat)
  rbeta(1, shape1 = s + 1, shape2 = n - s + 1)
}

latent_f <- function(theta) {
  as.matrix(rbinom(4000, 1, theta), ncol = 1)
}

st_f <- function(i, xi, sdp) {
  x    <- rep(0, 4000)
  x[i] <- xi
  x
}

priv_f <- function(sdp, lx) {
  t1 <- sum(sdp == lx)
  t2 <- sum(sdp != lx)
  t1 * log(3/4) + t2 * log(1/4)
}


dmod <- new_privacy(post_f = post_f,
                    latent_f = latent_f,
                    priv_f = priv_f,
                    st_f = st_f,
                    npar = 1)


set.seed(1)
x_conf <- rbinom(4000, 1, .8)
idp <- as.logical(rbinom(4000, 1, .5))

x_sdp <- x_conf
x_sdp[idp] <- rbinom(sum(idp), 1, .5)

library(progressr)
handlers("progress")
with_progress({
out <- dapper_sample(dmod,
                     sdp = x_sdp,
                     warmup = 0,
                     init_par = .1,
                     niter = 400)
})

profvis({
dapper_chain(dmod,
              sdp = x_sdp,
              init_par = .8,
              niter = 10)

})

dmat <- latent_f(.25)
st_f(1, dmat)
reduce(lapply(1:nrow(dmat), function(i) st_f(i, dmat)), "+")




set.seed(1)
tmp <- apply(UCBAdmissions, 3, identity, simplify=FALSE)
adm_cnf <- Reduce('+', tmp)

tdf <- tibble(gender = c("Male", "Male", "Female", "Female"),
              status = c("Admitted", "Rejected", "Admitted", "Rejected"),
              n = c(1198, 1493, 557, 1278))


tdf <- tibble(gender = c(1,1,0,0),
              status = c(1,0,1,0),
              n = c(1198, 1493, 557, 1278))

set.seed(1)
ix1 <- as.logical(rbinom(4526, 1, .5))
tx1 <- sum(ix1)
rx1 <- rbinom(tx1, 1, .5)

ix2 <- as.logical(rbinom(4526, 1, .5))
tx2 <- sum(ix2)
rx2 <- rbinom(tx2, 1, .5)

tdf2 <- uncount(tdf, n)
tdf2$pgender <- tdf2$gender
tdf2$pgender[ix1] <- rx1

tdf2$pstatus <- tdf2$status
tdf2$pstatus[ix2] <- rx2

sdp <- c(tdf2$pgender, tdf2$pstatus)


tdf <- c(rep(c(1,0,0,0), 1198),
         rep(c(0,1,0,0), 1493),
         rep(c(0,0,1,0), 557),
         rep(c(0,0,0,1), 1278))

m1 <- matrix(0, nrow = 1198, ncol = 4)
m1[,1] <- 1

m2 <- matrix(0, nrow = 1493, ncol = 4)
m2[,2] <- 1

m3 <- matrix(0, nrow = 557, ncol = 4)
m3[,3] <- 1

m4 <- matrix(0, nrow = 1278, ncol = 4)
m4[,4] <- 1

m <- rbind(m1,m2,m3,m4)
msdp <- m
ix <- as.logical(rbinom(nrow(msdp), 1, .50))
for(i in 1:length(ix)) {
  if(ix[i]){
    msdp[i,] <- 0
    r <- sample(1:4,1)
    msdp[i,r] <- 1
  }
}

sdp <- c(msdp[,1], msdp[,2], msdp[,3], msdp[,4])
#----------------

post_f <- function(dmat, theta) {
  x <- apply(dmat,2,sum)
  t1 <- rgamma(4, x + 1, 1)
  t1/sum(t1)
}


latent_f <- function(theta) {
 t(rmultinom(4526, 1, theta))
}

st_f <- function(i, xi, sdp) {
  x    <- rep(0, 4526 * 4)
  x[i] <- xi[1]
  x[i + 4526] <- xi[2]
  x[i + 2 * 4526] <- xi[3]
  x[i + 3 * 4526] <- xi[4]
  x
}

priv_f <- function(sdp, tx) {
  m1 <- matrix(sdp, nrow = 4526, ncol = 4, byrow = FALSE)
  m2 <- matrix(tx, nrow = 4526, ncol = 4, byrow = FALSE)
  m3 <- m1 == m2
  m4 <- m3[,1] & m3[,2] & m3[,3] & m3[,4]
  t1 <- sum(m4)
  t1 * log(5/8) + (4526 - t1) * log(3/8)
}

priv_f <- function(sdp, tx) {
  m1 <- matrix(sdp, nrow = 4526, ncol = 4, byrow = FALSE)
  m2 <- matrix(tx, nrow = 4526, ncol = 4, byrow = FALSE)
  m3 <- m1 == m2
  m4 <- m3[,1] & m3[,2] & m3[,3] & m3[,4]
  t1 <- sum(m4)
  t1 * log(5/8) + (4526 - t1) * log(3/8)
}


dmod <- new_privacy(post_f = post_f,
                    latent_f = latent_f,
                    priv_f = priv_f,
                    st_f = st_f,
                    npar = 4)


library(progressr)
with_progress({
  out <- dapper_sample(dmod,
                       sdp = x_sdp,
                       init_par = rep(.25,4),
                       niter = 5)
})

library(profvis)
profvis({
dapper_chain(dmod,
             sdp = sdp,
             init_par = rep(.25,4),
             niter = 5)
})


dapper_chain(dmod,
             sdp = sdp,
             init_par = rep(.25,4),
             niter = 10)

t1 <- apply(sapply(1:nrow(m), function(i) st_f(i, msdp[i,], sdp)), 1, sum)


#-------------------------------------------------------------------------------

post_f <- function(dmat, theta) {
  x <- c(table(as.data.frame(dmat)))
  t1 <- rgamma(4, x + 1, 1)
  t1/sum(t1)
}


latent_f <- function(theta) {
  tl <- list(c(1,1), c(1,0), c(0,1), c(0,0))
  rs <- sample(tl, 200, replace = TRUE, prob = theta)
  do.call(rbind, rs)
}

st_f <- function(i, xi, sdp) {
  x    <- rep(0, 200 * 2)
  x[i] <- xi[1]
  x[i + 200] <- xi[2]
  x
}

priv_f <- function(sdp, tx) {
  t1 <- sum(sdp == tx)
  t1 * log(3/4) + (400 - t1) * log(1/4)
}

set.seed(1)
tmp <- apply(UCBAdmissions, 3, identity, simplify=FALSE)
adm_cnf <- Reduce('+', tmp)

tdf <- tibble(gender = c(1,1,0,0),
              status = c(1,0,1,0),
              n = c(1198, 1493, 557, 1278)) %>% uncount(n)

tdf <- tdf[sample(1:nrow(tdf), 200),]


set.seed(1)
sdp <- c(tdf$gender, tdf$status)
ix <- as.logical(rbinom(200 * 2, 1, .5))
sdp[ix] <- rbinom(sum(ix), 1, .5)


dmod <- new_privacy(post_f = post_f,
                    latent_f = latent_f,
                    priv_f = priv_f,
                    st_f = st_f,
                    npar = 4)


library(progressr)
library(furrr)
handlers(global = TRUE)
plan(multisession, workers = 4)

handlers("cli")
out <- dapper_sample(dmod,
                       sdp = sdp,
                       warmup = 200,
                       chains = 4,
                       init_par = c(.25, .10, .15, .5),
                       niter = 2000)


library(profvis)
profvis({
  dapper_chain(dmod,
               sdp = sdp,
               init_par = rep(.25,4),
               niter = 5)
})



out <- dapper_chain(dmod,
             sdp = sdp,
             init_par = c(.25, .10, .15, .5),
             niter = 4)
