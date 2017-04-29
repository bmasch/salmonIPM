functions {
  # spawner-recruit functions
  real SR(real a, real b, real S, real A) {
    real R;
    R = a * S / (A + b * S);
    return(R);
  }
  
  # Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(real y, real mu, real sigma, real shape) {
    return(-(fabs(y - mu)/sigma)^shape);
  }
}

data {
  int<lower=1> N;               # total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];  # population identifier
  int<lower=1,upper=N> year[N]; # brood year identifier
  int<lower=1,upper=N> N_fit;   # number of cases used in fitting (non-missing S and R)
  int<lower=1,upper=N> which_fit[N_fit]; # cases used in fitting
  vector<lower=0>[N] S;         # observed annual total spawner abundance (not density)
  vector<lower=0>[N] R;         # total natural recruit abundance (not density), including harvest and broodstock removals
  vector[N] A;                  # habitat area associated with each spawner abundance obs
  int<lower=0,upper=1> S_NA[N]; # logical indicating whether S is missing and should be simulated
  int<lower=0,upper=1> R_NA[N]; # logical indicating whether R is missing and should be simulated
  int<lower=2> N_age;           # number of adult age classes
  int<lower=2> max_age;         # maximum adult age
  matrix<lower=0,upper=1>[max(pop),N_age] p;  # average recruit age distributions for each pop 
}

transformed data {
  int<lower=1,upper=N> N_pop;          # number of populations
  int<lower=1,upper=N> N_year;         # number of years
  int<lower=2> ages[N_age];              # adult ages
  
  N_pop = max(pop);
  N_year = max(year);
  for(j in 1:N_age)
    ages[j] = max_age - N_age + j;
}

parameters {
  vector<lower=0>[N_pop] a;             # intrinsic productivity
  vector<lower=0>[N_pop] b;             # density dependence
  vector<lower=-1,upper=1>[N_pop] rho;  # AR(1) coefs of residuals
  vector<lower=0>[N_pop] sigma;         # residual error SD
}

transformed parameters {
  vector<lower=0>[N] R_hat;             # expected recruit abundance (not density) by brood year
  vector<lower=0>[N] R_ar1;             # expected recruit abundance, taking AR(1) errors into account
  vector<lower=0>[N] sigma_ar1;         # residual error SD for each observation
  
  # Predict recruitment
  R_hat = rep_vector(0,N);
  R_ar1 = rep_vector(0,N);
  sigma_ar1 = rep_vector(0,N);
    
  for(i in 1:N_fit)
  {
    R_hat[which_fit[i]] = A[which_fit[i]] * SR(a[pop[which_fit[i]]], b[pop[which_fit[i]]], S[which_fit[i]], A[which_fit[i]]);
    if(i==1 || pop[which_fit[i-1]] != pop[which_fit[i]])
    {
      R_ar1[which_fit[i]] = R_hat[which_fit[i]];
      sigma_ar1[which_fit[i]] = sigma[pop[which_fit[i]]];
    }
    else
    {
      real err;    # temp variable: residual at the last non-missing observation
      int dt;      # temp variable: number of years since last non-missing observation
      real rho2j;  # temp variable: sum of powers of rho
      
      err = log(R[which_fit[i-1]]) - log(R_hat[which_fit[i-1]]);
      dt = which_fit[i] - which_fit[i-1];
      R_ar1[which_fit[i]] = R_hat[which_fit[i]]*exp(rho[pop[which_fit[i]]]^dt * err);
      
      rho2j = 0;
      for(j in 0:(dt-1))
        rho2j = rho2j + rho[pop[which_fit[i]]] ^ (2*j);
      sigma_ar1[which_fit[i]] = sigma[pop[which_fit[i]]] * sqrt(rho2j);
    }
  }
}

model {
  # Priors
  a ~ lognormal(0,5);
  b ~ lognormal(0,10);
  for(i in 1:N_pop)
  {
    rho[i] ~ pexp(0,0.8,10);   # mildly regularize rho to ensure stationarity
    sigma[i] ~ pexp(0,2,10);
  }

  # Likelihood
  R[which_fit] ~ lognormal(log(R_ar1[which_fit]), sigma_ar1[which_fit]);
}

generated quantities {
  vector[N] S_sim;    # simulated spawners
  vector[N] R_sim;    # simulated recruits
  vector[N] err_sim;  # simulated AR(1) residual errors
  
  S_sim = S;
  R_sim = R;
  err_sim = rep_vector(0,N);
  
  for(i in 1:N)
  {
    if(S_NA[i] == 1)
    {
      if(i >= max_age && pop[i-max_age] == pop[i])
      {
        S_sim[i] = 0;
        for(j in 1:N_age)
          S_sim[i] = S_sim[i] + R_sim[i-ages[j]]*p[pop[i],j];
      }
    }
    
    if(i == 1 || pop[i-1] != pop[i])
      err_sim[i] = normal_rng(0, sigma[pop[i]]);
    else
      err_sim[i] = normal_rng(rho[pop[i]]*err_sim[i-1], sigma[pop[i]]);
    
    if(R_NA[i] == 1)
      R_sim[i] = A[i]*SR(a[pop[i]], b[pop[i]], S_sim[i], A[i])*exp(err_sim[i]);
  }
}
