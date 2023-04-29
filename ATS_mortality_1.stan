//  ATS_mortality.stan
// Scott Akenhead scott@s4s.com 2022-12-03

data {
  int<lower=0>  N;               // number observations, 116
  int<lower=0>  n_years;         // number years, 24
  int<lower=0>  n_obs[n_years];  // number obs each year
  vector [N]    day ;            // day of each obs, -151 to 188
  vector<lower=0> [N] y ;        // abundance, ATS estimate
  vector<lower=0> [N] prec;      // 0 > regression_weight <= 1
}
parameters {
  real<lower=0> m_0;               // mortality rate, per day.
  vector <lower=0> [n_years] N_0;  // fitted abundance day
  real<lower=0> sigma;             // std deviation of residuals
}
transformed parameters{
  real m0;
  m0 = m_0 *0.001;                 // sample about 5, use as 0.005
}
model {
  int n1; int n2;
  vector [N] mt;
  vector [N] y_hat;                      // predictions of observations
  mt = exp(-(m0 * day));                 // survival to each obs day in all years
  n1 = 1;
  for (j in 1:n_years) {                 // each smolt_year
    n2    = n1 + n_obs[j] -1 ;
    y_hat[n1:n2] = N_0[j] * mt[n1:n2];  // predictions
    n1    = n2+1;                       // for next obs
  }
  target += prec * normal_lpdf(y | y_hat, sigma); // log likelihood
}
