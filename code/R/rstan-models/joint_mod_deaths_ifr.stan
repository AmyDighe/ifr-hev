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
  
  // Prob HEV death given jaundice death
  vector[L] alpha_d = alpha_mu_d + alpha_sigma_d * alpha_std_d;    // place random intercepts
  vector[J] p_HEV_d = inv_logit(alpha_d[map_to_place] + X*beta_d);    // probability of being HEV positive given fatal jaundice
  
    
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
  alpha_mu ~ normal(-6, 2);          // alpha_mu ~ normal(-7.6, 1);  mean IFR ~ 0.05%
  alpha_sigma ~ normal(0, 0.5);
  beta ~ normal(0, 0.5);
  
  // Priors for prob of HEV death given jaundice death
  alpha_std_d ~ std_normal();
  alpha_mu_d ~ normal(-2, 1);
  alpha_sigma_d ~ normal(0, .5);
  beta_d ~ normal(0, 0.5);
  
  // Latent infections
  I_std ~ std_normal(); // automatically vectorized over j
  
  // Loglik of HEV deaths among jaundice deaths in hospitals
  DHEV ~ binomial(DAJ, p_HEV_d);
  
  // Likelihood: propagate uncertainty in deaths
  for (j in 1:J) {
    int n_max_d = DcAJ[j];    // maximum number of HEV deaths (= all jaundice deaths in community)
    vector[n_max_d] lp;
    real prob;
    prob = fmax((p_HEV_d[j] + spec - 1) / (spec + sens - 1), 0);    // Equivalent to: if (prob >= 0) then use prob, else 0
    prob = fmin(prob, 1);    // Equivalent to: if (prob < 1) then use prob, else 1
    
    for (d in 1:n_max_d) {
      // Poisson likelihood scaled by latent infections
      lp[d] = binomial_lpmf(d | n_max_d, prob) + poisson_lpmf(d | I_latent[j] * p[j]);
    }
    target += log_sum_exp(lp);
  }
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
