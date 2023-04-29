//  ATS_mortality_2.stan
// Scott Akenhead scott@s4s.com 2022-12-04
functions{
  vector surv (real m_0, real m_1, int n, vector all_days, int [] day_index) {
  vector [n] surv;
  surv = exp(-m_0 * exp(-m_1 * all_days)) ;
  for(j in 2:n) surv[j] = surv[j] * surv[j-1]; // cumulative product
  surv = surv / surv[151];    // survival at day 0 is 1
  return surv[day_index];     // survival from day 0 to each obs day.
  }
}
data {
  int<lower=0>        N;               // n observations, 116
  int<lower=0>        n_years;         // n  years, 24
  int<lower=0>        n_obs[n_years];  // n obs each year
  vector [N]          day ;            // day of each obs, -151 to 188
  vector<lower=0> [N] y ;              // abundance, ATS obs
  vector<lower=0> [N] prec;            // regression weights
                                       // for function surv()
  vector [340]        all_days ;       // sequence -151 to 188
  int                 day_index[N] ;   // index obs days in all_days
}
parameters {
  real<lower=0> m0;                   // mortality rate
  real<lower=0> m1;                  // rate of decline of m_0
  vector <lower=0> [n_years] N_0;     // fitted abundance day
  real<lower=0> sigma;                // std deviation of residuals
}
transformed parameters{
  real<lower=0> m_0;
  real<lower=0> m_1;
  m_0 = m0 * 0.001;
  m_1 = m1 * 0.001;
}
model {
  int n1; int n2;
  vector [N] y_hat;                   // predictions of observations
  vector [N] survive;                 // survival day zero to day of obs
  survive = surv(m_0, m_1, 340, all_days, day_index);
  n1 = 1;
  for (j in 1:n_years) {              // each smolt_year
    n2    = n1 + n_obs[j] -1 ;
    y_hat[n1:n2] = N_0[j] * survive[n1:n2];  // predictions
    n1    = n2+1;                     // for next obs
  }
  target += prec * normal_lpdf(y | y_hat, sigma);
}
