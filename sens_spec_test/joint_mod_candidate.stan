// CANDIDATE replacement for joint_mod_deaths_ifr.stan (exploratory).
//  1. Test error applied forwards to the hospital data (no clipping).
//  2. Community jaundice deaths are the observed data:
//       DcAJ ~ Poisson(I * IFR / p_HEV_d)
//     in place of the sum over d of Binomial(d | DcAJ, p_HEV_d) * Poisson(d | I * IFR).
//  3. All priors are passed as data so they can be varied without recompiling.
//  4. The HEV fraction has its own covariate matrix X_d, so a covariate (e.g. pregnancy)
//     can be dropped from it by zeroing its column.
data {
  int<lower=0> J;
  int<lower=1> L;
  int<lower=0> M;
  array[J] int<lower=0> DHEV;   // HEV IgM+ deaths in hospital
  array[J] int<lower=0> DAJ;    // acute jaundice deaths tested in hospital
  array[J] int<lower=0> DcAJ;   // acute jaundice deaths recalled by community
  vector[J] I_post_mean;
  vector[J] I_post_sd;
  matrix[J, M] X;               // covariates for IFR
  matrix[J, M] X_d;             // covariates for the HEV fraction
  array[J] int<lower=1, upper=L> map_to_place;
  real<lower=0, upper=1> sens;
  real<lower=0, upper=1> spec;
  // priors
  real alpha_mu_mean;  real<lower=0> alpha_mu_sd;  real<lower=0> alpha_sigma_sd;
  real alpha_mu_d_mean; real<lower=0> alpha_mu_d_sd; real<lower=0> alpha_sigma_d_sd;
  vector[M] beta_mean;   vector<lower=0>[M] beta_sd;
  vector[M] beta_d_mean; vector<lower=0>[M] beta_d_sd;
}
parameters {
  vector[L] alpha_std;
  real alpha_mu;
  real<lower=1e-6> alpha_sigma;
  vector[M] beta;
  vector[L] alpha_std_d;
  real alpha_mu_d;
  real<lower=1e-6> alpha_sigma_d;
  vector[M] beta_d;
  vector[J] I_std;
}
transformed parameters {
  vector[L] alpha = alpha_mu + alpha_sigma * alpha_std;
  vector[J] p = inv_logit(alpha[map_to_place] + X * beta);              // IFR
  vector[L] alpha_d = alpha_mu_d + alpha_sigma_d * alpha_std_d;
  vector[J] p_HEV_d = inv_logit(alpha_d[map_to_place] + X_d * beta_d);  // TRUE HEV fraction
  vector[J] p_HEV_d_obs = sens * p_HEV_d + (1 - spec) * (1 - p_HEV_d);  // apparent IgM+ fraction
  vector[J] I_latent = I_post_mean .* exp((I_post_sd ./ I_post_mean) .* I_std);
}
model {
  alpha_std ~ std_normal();
  alpha_mu ~ normal(alpha_mu_mean, alpha_mu_sd);
  alpha_sigma ~ normal(0, alpha_sigma_sd);
  beta ~ normal(beta_mean, beta_sd);
  alpha_std_d ~ std_normal();
  alpha_mu_d ~ normal(alpha_mu_d_mean, alpha_mu_d_sd);
  alpha_sigma_d ~ normal(0, alpha_sigma_d_sd);
  beta_d ~ normal(beta_d_mean, beta_d_sd);
  I_std ~ std_normal();
  DHEV ~ binomial(DAJ, p_HEV_d_obs);
  DcAJ ~ poisson(I_latent .* p ./ p_HEV_d);
}
generated quantities {
  vector[J] deaths_exp_ifr = I_latent .* p;          // expected HEV deaths from infections x IFR
  vector[J] deaths_exp_phevd;                        // observed jaundice deaths x HEV fraction
  array[J] int DcAJ_rep = poisson_rng(I_latent .* p ./ p_HEV_d);
  array[J] int DHEV_rep = binomial_rng(DAJ, p_HEV_d_obs);
  for (j in 1:J) deaths_exp_phevd[j] = DcAJ[j] * p_HEV_d[j];
}
