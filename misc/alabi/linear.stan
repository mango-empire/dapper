data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;

}

parameters {
  real<lower=0> sigma;
  real beta;
  real alpha;
}

model {
  for(i in 1:N) {
    y[i] ~ normal(alpha + x[i] * beta, sigma);
  }
}
