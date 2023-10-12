library(tidyverse)
library(stringr)

fnames <- list.files("data_out\\data")

hrd <- function(fn) {
  tmp <- str_split(fn, "_")[[1]]
  eps <- as.numeric(tmp[2])
  run <- as.integer(tmp[3])
  tdf <- read_csv(paste0("data_out\\data\\",fn))
  tdf$epsilon <- eps
  tdf$run <- run
  tdf
}

tmp <- map(fnames, function(fn) hrd(fn))

small_df <- do.call(rbind, tmp)
write_csv(small_df, "small_df.csv")

big_df$niter <- 95000
small_df$niter <- 5000

full_df <- rbind(big_df,small_df)
write_csv(full_df, "full_df.csv")
full_df <- read_csv("full_df.csv")



full_df %>%
  mutate(epsilon = as_factor(epsilon), niter = as_factor(niter)) %>%
  group_by(variable, epsilon) %>%
  ggplot(aes(x=epsilon, y=mean, color=niter)) + geom_boxplot() + geom_jitter() + facet_wrap(~variable)

full_df %>%
  mutate(epsilon = as_factor(epsilon), niter = as_factor(niter)) %>%
  group_by(variable, epsilon) %>%
  ggplot(aes(x=epsilon, y=rhat, color=niter)) + geom_boxplot() + facet_wrap(~variable)



full_df %>%
  mutate(epsilon = as_factor(epsilon), niter = as_factor(niter)) %>%
  group_by(variable, epsilon) %>%
  ggplot(aes(x=epsilon, y=ess_bulk, color=niter)) + geom_boxplot() + facet_wrap(~variable)


big_df %>%
  mutate(epsilon = as_factor(epsilon)) %>%
  group_by(variable, epsilon) %>%
  ggplot(aes(epsilon, ess_bulk, color = variable)) + geom_boxplot() + geom_jitter() + facet_wrap(~variable)



fnames <- list.files()

hrd <- function(fn) {
  tmp <- str_split(fn, "_")[[1]]
  eps <- as.numeric(tmp[2])
  run <- as.integer(tmp[3])
  tdf <- read_csv(fn) %>% select(variable, mean)
  tdf$epsilon <- eps
  tdf$run <- run
  tdf
}

tmp <- map(fnames, function(fn) hrd(fn))

beta_df <- do.call(rbind, tmp)

write_csv(beta_df, "beta_df.csv")

beta_df %>%
  mutate(epsilon = as_factor(epsilon)) %>%
  group_by(variable, epsilon) %>%
  ggplot(aes(epsilon, mean, color = variable)) + geom_boxplot() + geom_jitter() + facet_wrap(~variable)

#----------


hrd <- function(fn) {
  tmp <- str_split(fn, "_")[[1]]
  eps <- as.numeric(tmp[2])
  run <- as.integer(tmp[3])
  tdf <- read_csv(fn) %>% select(variable, rhat)
  tdf$epsilon <- eps
  tdf$run <- run
  tdf
}

tmp <- map(fnames, function(fn) hrd(fn))

rhat_df <- do.call(rbind, tmp)

write_csv(rhat_df, "rhat_df.csv")

rhat_df %>%
  mutate(epsilon = as_factor(epsilon)) %>%
  group_by(variable, epsilon) %>%
  ggplot(aes(epsilon, rhat, color = variable)) + geom_boxplot() + geom_jitter() + facet_wrap(~variable)


#------------

hrd <- function(fn) {
  tmp <- str_split(fn, "_")[[1]]
  eps <- as.numeric(tmp[2])
  run <- as.integer(tmp[3])
  tdf <- read_csv(fn) %>% select(variable, ess_bulk)
  tdf$epsilon <- eps
  tdf$run <- run
  tdf
}

tmp <- map(fnames, function(fn) hrd(fn))

ess_bulk_df <- do.call(rbind, tmp)

write_csv(ess_bulk_df, "ess_bulk_df.csv")

ess_bulk_df %>%
  mutate(epsilon = as_factor(epsilon)) %>%
  group_by(variable, epsilon) %>%
  ggplot(aes(epsilon, ess_bulk, color = variable)) + geom_boxplot() + geom_jitter() + facet_wrap(~variable)
