functions {
  
  // Continuous approximation to the Negative Binomial (Gamma–Poisson mixture)
  // Replace the binomial coefficient with ratio of gamma functions
  // Following: Javier and https://stats.stackexchange.com/questions/310676/continuous-generalization-of-the-negative-binomial-distribution
  real neg_binomial_2_cont_lpdf(real x, real mu, real phi) {
    real log_gamma_ratio;
    real ll;

    log_gamma_ratio = lgamma(x + phi) - lgamma(x + 1) - lgamma(phi);
    ll = log_gamma_ratio + x * log(mu) - phi * log1p(mu / phi) - x * log(mu + phi);

    return ll;
  }
  
  
  // point estimates for mid-point of age class
  
    // model 1 - constant FoI
    
    real p_seroprev_1_point(real age,
                            real foi){
            real pred_seroprev;
            
            pred_seroprev = 1 - exp(-foi * age);
            
            return pred_seroprev;
                            }
  
    // model 1 - age-varying FoI
    
    real p_seroprev_1varFoI_point(real age,
                                  real age_foi,
                                  real foi1,
                                  real foi2){
            real pred_seroprev;
            
            if(age <= age_foi){
              pred_seroprev = 1 - exp(-foi1 * age);
                            } else {
              pred_seroprev = 1 - exp(-foi1 * age_foi - foi2 * (age - age_foi));
              }
            
              
              return pred_seroprev;
            }
            
  
    // model 2 - constant FoI, constant sigma
    
    real p_seroprev_2_point(real age,
                            real foi,
                            real sigma){
            real pred_seroprev;
            
            pred_seroprev = (foi / (foi + sigma)) * (1 - exp(-(foi + sigma) * age));
            
            return pred_seroprev;
                            }
  
    // model 2 - age-varying FoI, constant sigma
    
    real p_seroprev_2_varFoI_point(real age,
                                   real age_foi,
                                   real foi1,
                                   real foi2,
                                   real sigma){
            real pred_seroprev;
            
            if(age <= age_foi){
              pred_seroprev = p_seroprev_2_point(age, foi1, sigma);
            } else{
              pred_seroprev = ((foi1 / (foi1 + sigma)) * (1 - exp(-(foi1 + sigma) * age)) - (foi2 / (foi2 + sigma))) *
              exp(-(foi2 + sigma) * (age - age_foi)) + (foi2 / (foi2 + sigma));
            }
            
            return pred_seroprev;
                                   }
  
    // model 2 - age varying FoI, age-varying sigma - when age_sigma < age < age_foi
  
    real p_seroprev_2ii_point(real age,
                              real age_sigma,
                              real foi1,
                              real sigma1,
                              real sigma2){
        real pred_seroprev;
        
        pred_seroprev = 1 - ((1 - p_seroprev_2_point(age_sigma, foi1, sigma1)) *
        (1 - (foi1/(foi1 + sigma2)) * (1 - exp(-(foi1 + sigma2) * (age - age_sigma)))));
        
        return pred_seroprev;
                          }
    // model 2 - age varying FoI, age-varying sigma - when age > age_foi
    
    real p_seroprev_2iii_point(real age,
                               real age_sigma,
                               real age_foi,
                               real foi1,
                               real foi2,
                               real sigma1,
                               real sigma2){
        real pred_seroprev;
        
        pred_seroprev = 1 - ((1 - p_seroprev_2ii_point(age_foi, age_sigma, foi1, sigma1, sigma2)) *
        (1 - (foi2 / (foi2 + sigma2)) * (1 - exp(-(foi2 + sigma2)*(age - age_foi)))));
        
        return pred_seroprev;
                               }
                          
  // estimates for average across the age class spanning age1 - age2

  real p_seroprev_1(real age1,
                    real age2, 
                    real foi){
                      
        real pred_seroprev;
        
        pred_seroprev = (1 / (age2 - age1)) * (age2 - age1 + (1 / foi) * (exp(-foi * age2) - exp(-foi * age1)));
        
        return pred_seroprev;
    
  }
  
  real p_seroprev_1_varFoI(real age1,
                          real age2,
                          real age_foi,
                          real foi1,
                          real foi2){
      
      real pred_seroprev;
                            
      pred_seroprev = (1 / (age2 - age1)) * 
      (age2 - age1 + (1 / foi2) * (exp(-(foi1 * age_foi + foi2 * (age2 - age_foi))) - exp(-(foi1 * age_foi + foi2 * (age1 - age_foi)))));
    
                            
      return pred_seroprev;
  }
  
  real p_seroprev_2(real age1,
                    real age2, 
                    real foi,
                    real sigma){
                      
       real pred_seroprev;
       
       pred_seroprev = (1 / (age2 - age1)) * 
                        (foi / (foi + sigma)) * 
                        ((age2 - age1) + (1 / (foi + sigma)) * (exp(-(foi + sigma) * age2) - exp(-(foi + sigma) * age1)));
      
      return pred_seroprev;
  }
  
  real p_seroprev_2_varFoI(real age1,
                          real age2,
                          real age_foi,
                          real foi1,
                          real foi2,
                          real sigma){
                            
      real pred_seroprev;
                            
      pred_seroprev = (1 / (age2 - age1)) * ((age2 - age1) * (foi2 / (foi2 + sigma)) 
                      + (1 / (foi2 + sigma)) * ((foi1/(foi1 + sigma)) * (1 - exp(-(foi1 + sigma) * age_foi)) - (foi2/ (foi2 + sigma))) *
                      (exp(-(foi2 + sigma) * (age1 - age_foi)) - exp(-(foi2 + sigma) * (age2 - age_foi))));
                      
      return pred_seroprev;
      
  }
  
    real p_seroprev_2_varSigma(real age1,
                          real age2,
                          real age_sigma,
                          real foi1,
                          real sigma1,
                          real sigma2){
                            
      real pred_seroprev;
                            
      pred_seroprev = (1 / (age2 - age1)) *
      ((age2 - age1) - (1 - p_seroprev_2_point(age_sigma, foi1, sigma1)) * 
      ((age2 - age1) * (1 - (foi1 / (foi1 + sigma2))) - 
      (1 / (foi1 + sigma2)) * (foi1 / (foi1 + sigma2)) * exp((foi1 + sigma2) * age_sigma) *
      (exp(-(foi1 + sigma2) * age2) - exp(-(foi1 + sigma2) * age1))));
                      
      return pred_seroprev;
      
  }
  
      real p_seroprev_3_varSigma(real age1,
                          real age2,
                          real age_foi,
                          real age_sigma,
                          real foi1,
                          real foi2,
                          real sigma1,
                          real sigma2){
                            
      real pred_seroprev;
                            
      pred_seroprev = (1 / (age2 - age1)) *
      ((age2 - age1) - (1 - p_seroprev_2ii_point(age_foi, age_sigma, foi1, sigma1, sigma2)) * 
      ((age2 - age1) * (1 - (foi2 / (foi2 + sigma2))) - 
      (1 / (foi2 + sigma2)) * (foi2 / (foi2 + sigma2)) * exp((foi2 + sigma2) * age_foi) *
      (exp(-(foi2 + sigma2) * age2) - exp(-(foi2 + sigma2) * age1))));
                      
      return pred_seroprev;
      
  }
  
  // likelihoods - integrating across age class to get mean sp
  
  real cat1_lpmf(int seropos,
                int N,
                real foi,
                real age1,
                real age2){
    real pred_seroprev;
    real loglik;
    
    pred_seroprev = p_seroprev_1(age1, age2, foi);
    
    loglik = binomial_lpmf(seropos|N, pred_seroprev);
    
    return loglik;
  }

real cat2_lpmf(int seropos,
               int N,
               real foi,
               real sigma,
               real age1,
               real age2){
  real pred_seroprev;
  real loglik;
  
  pred_seroprev = p_seroprev_2(age1, age2, foi, sigma);
  
  loglik = binomial_lpmf(seropos|N, pred_seroprev);
  
  return loglik;
}

  real cat1_varFoI_lpmf(int seropos,
                int N,
                real foi1,
                real foi2,
                real age_foi,
                real age1,
                real age2){
    real pred_seroprev;
    real loglik;
    
    if(age2 <= age_foi){ //ensure no age classes span ages both above and beyond age_foi or there will be trouble
    pred_seroprev = p_seroprev_1(age1, age2, foi1);
    
  } else {
      pred_seroprev = p_seroprev_1_varFoI(age1, age2, age_foi, foi1, foi2);
    }
    

    loglik = binomial_lpmf(seropos|N, pred_seroprev);
    
    return loglik;
  }
  
    real cat2_varFoI_lpmf(int seropos,
                int N,
                real foi1,
                real foi2,
                real sigma,
                real age_foi,
                real age1,
                real age2){
    real pred_seroprev;
    real loglik;
    
    if(age2<=age_foi){//ensure no age classes span ages both above and beyond age_foi or there will be trouble
      
      pred_seroprev = p_seroprev_2(age1, age2, foi1, sigma);
      
    } else{
          pred_seroprev = p_seroprev_2_varFoI(age1, age2, age_foi, foi1, foi2, sigma);
          }
    
    loglik = binomial_lpmf(seropos|N, pred_seroprev);
    
    return loglik;
  }
  
      real cat2_varFoIvarSigma_lpmf(int seropos,
                int N,
                real foi1,
                real foi2,
                real sigma1,
                real sigma2,
                real age_foi,
                real age_sigma,
                real age1,
                real age2){
    real pred_seroprev;
    real loglik;
    
    if(age2<=age_sigma){
      
      pred_seroprev = p_seroprev_2(age1, age2, foi1, sigma1);
      
    } else if(age2>age_sigma && age2<=age_foi){
      
          pred_seroprev = p_seroprev_2_varSigma(age1, age2, age_sigma, foi1, sigma1, sigma2);
          
          } else if (age2>age_foi){
            
            pred_seroprev = p_seroprev_3_varSigma(age1, age2, age_foi, age_sigma, foi1, foi2, sigma1, sigma2);
            
          } 
    
    loglik = binomial_lpmf(seropos|N, pred_seroprev);
    return loglik;
  }
  
  // likelihoods - estimating sp at the mid-point of the age class
  
    real cat1point_lpmf(int seropos,
                int N,
                real foi,
                real age){
    real pred_seroprev;
    real loglik;
    
    pred_seroprev = p_seroprev_1_point(age, foi);
    
    loglik = binomial_lpmf(seropos|N, pred_seroprev);
    
    return loglik;
  }

real cat2point_lpmf(int seropos,
               int N,
               real foi,
               real sigma,
               real age){
  real pred_seroprev;
  real loglik;
  
  pred_seroprev = p_seroprev_2_point(age, foi, sigma);
  
  loglik = binomial_lpmf(seropos|N, pred_seroprev);
  
  return loglik;
}

  real cat1_varFoIpoint_lpmf(int seropos,
                int N,
                real foi1,
                real foi2,
                real age_foi,
                real age){
    real pred_seroprev;
    real loglik;
    
    if(age <= age_foi){
    pred_seroprev = p_seroprev_1_point(age, foi1);
    
  } else {
      pred_seroprev = p_seroprev_1varFoI_point(age, age_foi, foi1, foi2);
    }
    

    loglik = binomial_lpmf(seropos|N, pred_seroprev);
    
    return loglik;
  }
  
    real cat2_varFoIpoint_lpmf(int seropos,
                int N,
                real foi1,
                real foi2,
                real sigma,
                real age_foi,
                real age){
    real pred_seroprev;
    real loglik;
    
    if(age<=age_foi){
      
      pred_seroprev = p_seroprev_2_point(age, foi1, sigma);
      
    } else{
          pred_seroprev = p_seroprev_2_varFoI_point(age, age_foi, foi1, foi2, sigma);
          }
    
    loglik = binomial_lpmf(seropos|N, pred_seroprev);
    
    return loglik;
  }
  
    real cat2_varFoIvarSigmapoint_lpmf(int seropos,
                int N,
                real foi1,
                real foi2,
                real sigma1,
                real sigma2,
                real age_foi,
                real age_sigma,
                real age){
    real pred_seroprev;
    real loglik;
    
    if(age<=age_sigma){
      
      pred_seroprev = p_seroprev_2_point(age, foi1, sigma1);
      
    } else if(age>age_sigma && age<=age_foi){
      
          pred_seroprev = p_seroprev_2ii_point(age, age_sigma, foi1, sigma1, sigma2);
          
          } else if (age>age_foi){
            
            pred_seroprev = 1 - ((1 - p_seroprev_2ii_point(age_foi, age_sigma, foi1, sigma1, sigma2)) *
            (1 - (foi2/(foi2 + sigma2)) * (1 - exp(-(foi2+sigma2)*(age - age_foi)))));
            
          } 
    
    loglik = binomial_lpmf(seropos|N, pred_seroprev);
    return loglik;
  }
  
}
