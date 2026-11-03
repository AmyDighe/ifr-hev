#include functions.stan
data{
  int A; // number of age classes
  int Z; // number of posterior draws of seroprevalence from the INLA output
  real sp[A, Z]; // seroprevalence (posterior from INLA output)
  int age_foi; // cut off at which point FoI changes
  int age_sigma; // cut off at which point seroreversion rate changes
  vector[A] age1; // lower bound on age class (yrs)
  vector[A] age2; // upper bound on age class (yrs)
  real sigma1; // fix seroreversion rate for those < age_sigma based on longitudinal data
  real sigma2; // fix seroreversion rate for those > age_sigma based on longitudinal data
}

parameters{
  real <lower = 0> foi1; // Force of Infection
  real <lower = 0> foi2; // Force of Infection
  real <lower = 0> sd; // standard deviation of the logit-logit normal 
}

model{
  sd ~ normal(0.1, 0.1);
  foi1 ~ normal(0, 0.075);
  foi2 ~ normal(0, 0.075);
  real lp[A,Z];
  real pred_seroprev[A];

     for(a in 1:A){
       if(age2[a]<=age_sigma){
         pred_seroprev[a] = p_seroprev_2(age1[a], age2[a], foi1, sigma1);
         } else if(age2[a]>age_sigma && age2[a]<=age_foi){
           pred_seroprev[a] = p_seroprev_2_varSigma(age1[a], age2[a], age_sigma, foi1, sigma1, sigma2);
           } else if (age2[a]>age_foi){
             pred_seroprev[a] = p_seroprev_3_varSigma(age1[a], age2[a], age_foi, age_sigma, foi1, foi2, sigma1, sigma2);
             }
             for(z in 1:Z){
                 lp[a,z] = normal_lpdf(logit(sp[a,z])|logit(pred_seroprev[a]), sd);
                 }
     }

for(a in 1:A){
  target += log_sum_exp(to_vector(lp[a, ])) - log(Z);
}
  }


generated quantities{
  vector[A] seroprevalence;

  for(a in 1:A){
    
   if(age2[a] <= age_sigma){
      
      seroprevalence[a] = p_seroprev_2(age1[a], age2[a], foi1, sigma1);
      
    } else if(age2[a] > age_sigma && age2[a] <= age_foi){
      
          seroprevalence[a] = p_seroprev_2_varSigma(age1[a], age2[a], age_sigma, foi1, sigma1, sigma2);
          
          } else if (age2[a] > age_foi){
            
            seroprevalence[a] = p_seroprev_3_varSigma(age1[a], age2[a], age_foi, age_sigma, foi1, foi2, sigma1, sigma2);
            
          } 
  }
}


