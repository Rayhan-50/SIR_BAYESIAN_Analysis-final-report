functions {
  // SIR ODE system: y[1]=S, y[2]=I, y[3]=R
  vector sir_ode(real t,
                 vector y,
                 array[] real theta,
                 array[] real x_r,
                 array[] int x_i) {
    real beta  = theta[1];
    real gamma = theta[2];
    real N     = x_r[1];

    real S = y[1];
    real I = y[2];

    vector[3] dydt;
    dydt[1] = -beta * S * I / N;
    dydt[2] =  beta * S * I / N - gamma * I;
    dydt[3] =  gamma * I;
    return dydt;
  }
}

data {
  int<lower=1> N_weeks;
  array[N_weeks] int<lower=0> cases;
  real<lower=0> pop;

  // Prior hyperparameters — passed in as data for sensitivity analysis.
  // Baseline values (matching original model):
  //   R0_mu=log(2.0), R0_sigma=0.3
  //   gamma_mu=log(0.5), gamma_sigma=0.3
  //   rho_alpha=2, rho_beta=10
  //   phi_rate=0.5
  //   I0_mu=log(1), I0_sigma=2.0
  //   lambda_rate=0.5
  real          R0_mu;
  real<lower=0> R0_sigma;
  real          gamma_mu;
  real<lower=0> gamma_sigma;
  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;
  real<lower=0> phi_rate;
  real          I0_mu;
  real<lower=0> I0_sigma;
  real<lower=0> lambda_rate;
}

transformed data {
  array[1] real x_r = {pop};
  array[0] int  x_i;
  real t0 = 0.0;
  array[N_weeks] real ts;
  for (i in 1:N_weeks) ts[i] = i;
}

parameters {
  real<lower=0>          R0;
  real<lower=0>          gamma;
  real<lower=0, upper=1> rho;
  real<lower=0>          phi;
  real<lower=0.001>      I0;
  real<lower=0>          lambda;
}

transformed parameters {
  real beta = R0 * gamma;
  array[2] real theta = {beta, gamma};

  vector[3] y0;
  y0[1] = pop - I0;
  y0[2] = I0;
  y0[3] = 0.0;

  array[N_weeks] vector[3] y_sol;
  y_sol = ode_rk45(sir_ode, y0, t0, ts, theta, x_r, x_i);

  // Incidence = new infections per week = decrease in susceptibles (flow, not stock)
  array[N_weeks] real incidence;
  incidence[1] = fmax(1e-6, (pop - I0) - y_sol[1][1]);
  for (i in 2:N_weeks) {
    incidence[i] = fmax(1e-6, y_sol[i-1][1] - y_sol[i][1]);
  }

  // Expected cases: epidemic incidence scaled by reporting rate + background
  array[N_weeks] real mu;
  for (i in 1:N_weeks) {
    mu[i] = fmax(1e-6, rho * incidence[i] + lambda);
  }
}

model {
  R0     ~ lognormal(R0_mu,     R0_sigma);
  gamma  ~ lognormal(gamma_mu,  gamma_sigma);
  rho    ~ beta(rho_alpha,      rho_beta);
  phi    ~ exponential(phi_rate);
  I0     ~ lognormal(I0_mu,     I0_sigma);
  lambda ~ exponential(lambda_rate);

  for (i in 1:N_weeks) {
    cases[i] ~ neg_binomial_2(mu[i], phi);
  }
}

generated quantities {
  array[N_weeks] real cases_rep;
  array[N_weeks] real log_lik;

  real R0_prior     = lognormal_rng(log(2.0), 0.3);
  real gamma_prior  = lognormal_rng(log(0.5),  0.3);
  real rho_prior    = beta_rng(2, 10);
  real phi_prior    = exponential_rng(0.5);
  real I0_prior     = lognormal_rng(log(1), 2.0);
  real lambda_prior = exponential_rng(0.5);

  for (i in 1:N_weeks) {
    cases_rep[i] = neg_binomial_2_rng(mu[i], phi);
    log_lik[i]   = neg_binomial_2_lpmf(cases[i] | mu[i], phi);
  }
}
