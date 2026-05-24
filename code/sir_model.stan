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
  real<lower=0.001>      I0;      // initial seed — fractional allowed in ODE
  real<lower=0>          lambda;  // background sporadic cases per week
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

  // Incidence = new infections per week = decrease in susceptibles
  array[N_weeks] real incidence;
  incidence[1] = fmax(1e-6, (pop - I0) - y_sol[1][1]);
  for (i in 2:N_weeks) {
    incidence[i] = fmax(1e-6, y_sol[i-1][1] - y_sol[i][1]);
  }

  // Expected cases = epidemic incidence + background rate
  array[N_weeks] real mu;
  for (i in 1:N_weeks) {
    mu[i] = fmax(1e-6, rho * incidence[i] + lambda);
  }
}

model {
  // Priors
  R0     ~ lognormal(log(2.0),  0.3);   // median R0=2, allows 1.1–3.5
  gamma  ~ lognormal(log(0.5),  0.3);   // ~2-week infectious period
  rho    ~ beta(2, 10);                  // reporting ~17% mean
  phi    ~ exponential(0.5);
  I0     ~ lognormal(log(1),    2.0);   // wide prior: allows 0.001–100
  lambda ~ exponential(0.5);             // mean 2 background cases/wk

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
