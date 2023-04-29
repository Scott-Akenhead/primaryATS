//  Spawner Parr.stan
// Scott Akenhead scott@s4s.com 2022-11-26

data {
  int<lower=0>     n_years;  // number of years, cases
  vector[n_years]  spawner;  // x variable, predictor
  vector[n_years]  parr;     // y variable, to be predicted, predictee
}
parameters {
  real b_0;                 // intercept. Not regression through origin
  real b_1;                 // standard deviation in treatment effects
  real<lower=0> sigma;      // std deviation of residuals
}
model {
  parr ~ normal(b_0 + b_1 * spawner, sigma);
}
