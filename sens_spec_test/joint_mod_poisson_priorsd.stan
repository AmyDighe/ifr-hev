// PRIOR-SENSITIVITY variant: prior sd of covariate effects passed as data.
// EXPLORATORY variant (Poisson likelihood on DcAJ) of the joint model for P(HEV | fatal jaundice), HEV deaths in the community and IFR.
//
// Differs from joint_mod_deaths_ifr.stan in how imperfect test performance
// is handled: here the TRUE probability of HEV given fatal jaundice is the
// parameter, and test error is applied in the forward direction to get the
// apparent IgM+ probability that the hospital data inform. No back-correction
// (Rogan-Gladen) and therefore no clipping at 0/1 is needed.

data {
  int<lower=0> J;    // number of unique groups (n places x n age cats + n places)
  int<lower=1> L;          // number of places
  int<lower=0> M;          // number of covariates

  // All acute jaundice deaths and HEV IgM+ve deaths
  array[J] int<lower=0> DHEV;    // number of HEV IgM+ve deaths in hospital per group
  array[J] int<lower=0> DAJ;    // total number of acute jaundice deaths in hospital per group
  array[J] int<lower=0> DcAJ;   // total number of acute jaundice deaths recalled by community

  vector[J] I_post_mean;        // posterior mean of infections
  vector[J] I_post_sd;          // posterior sd of infections

  matrix[J, M] X;               // covariate matrix (no intercept)
  array[J] int<lower=1, upper=L> map_to_place; // group to place map

  // Imperfect test parameters
  real<lower=0, upper = 1> sens;
  real<lower=0, upper = 1> spec;

  // Smallest number of HEV deaths allowed among community jaundice deaths.
  // 1 reproduces the original model (sum over d = 1..DcAJ);
  // 0 also allows for no HEV deaths (sum over d = 0..DcAJ).
  int<lower=0, upper=1> d_min;

  // prior sd of the age/pregnancy effects (0.5 in the original model)
  real<lower=0> beta_sd;     // effects on logit IFR
  real<lower=0> beta_d_sd;   // effects on logit P(HEV | fatal jaundice)
}

parameters {
  vector[L] alpha_std;           // random intercepts (non-centered)
  real alpha_mu;                 // mean logit IFR
  real<lower=1e-6> alpha_sigma;     // sd of random intercepts
  vector[M] beta;                // covariate effects

  vector[L] alpha_std_d;    // place random intercept (non-centered) for prob of HEV death given jaundice death
  real alpha_mu_d;                // mean overall logit(p)
  real<lower=1e-6> alpha_sigma_d;    // sd of the overall p ~ logit normal
  vector[M] beta_d;    // covariate effects
  vector[J] I_std; // latent infections (standard normal reparameterisation)
}

transformed parameters {

  // IFR
  vector[L] alpha = alpha_mu + alpha_sigma * alpha_std;
  vector[J] p = inv_logit(alpha[map_to_place] + X * beta);      // IFR

  // TRUE prob of HEV given fatal jaundice
  vector[L] alpha_d = alpha_mu_d + alpha_sigma_d * alpha_std_d;    // place random intercepts
  vector[J] p_HEV_d = inv_logit(alpha_d[map_to_place] + X*beta_d);

  // APPARENT prob of testing IgM+ given fatal jaundice (true + false positives)
  vector[J] p_HEV_d_obs = sens * p_HEV_d + (1 - spec) * (1 - p_HEV_d);

  // Infections on the log-scale to avoid geometry issues
  vector[J] I_latent;
  for (j in 1:J) {
    I_latent[j] =
      I_post_mean[j] *
      exp( (I_post_sd[j] / I_post_mean[j]) * I_std[j] );
  }
}

model {

  // Priors
  alpha_std ~ std_normal();
  alpha_mu ~ normal(-6, 2);
  alpha_sigma ~ normal(0, 0.5);
  beta ~ normal(0, beta_sd);

  // Priors for (true) prob of HEV death given jaundice death
  alpha_std_d ~ std_normal();
  alpha_mu_d ~ normal(-2, 1);
  alpha_sigma_d ~ normal(0, .5);
  beta_d ~ normal(0, beta_d_sd);

  // Latent infections
  I_std ~ std_normal();

  // Loglik of HEV IgM+ deaths among jaundice deaths in hospitals:
  // the test result is what is observed, so use the apparent probability
  DHEV ~ binomial(DAJ, p_HEV_d_obs);

  // Community jaundice deaths as the observed data: if HEV deaths ~ Poisson(I * IFR)
  // and a fraction p_HEV_d of jaundice deaths are HEV, then all jaundice deaths
  // ~ Poisson(I * IFR / p_HEV_d)
  DcAJ ~ poisson(I_latent .* p ./ p_HEV_d);
}

generated quantities {
   // simulated deaths
  vector[J] deaths_sim_ifr;
  vector[J] deaths_sim_phevd;
  vector[J] deaths_exp_phevd;

  for (j in 1:J) {
    deaths_sim_ifr[j] = poisson_rng(I_latent[j] * p[j]);
    deaths_sim_phevd[j] = binomial_rng(DcAJ[j], p_HEV_d[j]);
    deaths_exp_phevd[j] = DcAJ[j]*p_HEV_d[j];
  }
}
