functions {
  # spawner-recruit functions
  real SR(real a, real Rmax, real S, real A) {
    real R;
    R = a*S/(A + a*S/Rmax);
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
  int<lower=1,upper=N> N_pop;      # number of populations
  int<lower=1,upper=N> N_year;     # number of years
  int<lower=2> ages[N_age];        # adult ages
  
  N_pop = max(pop);
  N_year = max(year);
  for(j in 1:N_age)
    ages[j] = max_age - N_age + j;
}

parameters {
  real mu_log_a;                        # hyper-mean log intrinsic productivity of wild spawners
  real<lower=0> sigma_log_a;            # hyper-SD log intrinsic productivity
  vector[N_pop] log_a_z;                # log intrinsic prod of wild spawners (Z-scores)
  real mu_log_Rmax;                     # hyper-mean log asymptotic recruitment
  real<lower=0> sigma_log_Rmax;         # hyper-SD log asymptotic recruitment
  vector[N_pop] log_Rmax_z;             # log asymptotic recruitment (Z-scores)
  real<lower=-1,upper=1> rho_log_aRmax; # correlation between log(a) and log(Rmax)
  real<lower=-1,upper=1> rho_log_phi;   # AR(1) coef for brood year log productivity anomalies
  real<lower=0> sigma_log_phi;          # hyper-SD of brood year log productivity anomalies
  vector[max(year)] log_phi_z;          # log brood year productivity anomalies (Z-scores)
  real<lower=0> sigma;                  # residual error SD
}

transformed parameters {
  vector<lower=0>[N_pop] a;     # intrinsic productivity 
  vector<lower=0>[N_pop] Rmax;  # asymptotic recruitment 
  vector<lower=0>[N_year] phi;  # brood year productivity anomalies
  vector<lower=0>[N] R_hat;     # expected recruit abundance (not density) by brood year
  
  # Multivariate Matt trick for [log(a), log(b)]
  {
    matrix[2,2] L_log_aRmax;      # Cholesky factor of correlation matrix of log(a) and log(Rmax)
    matrix[N_pop,2] aRmax;         # temp variable: matrix of a and Rmax
    vector[2] sigma_log_aRmax;     # temp variable: SD vector of [log(a), log(Rmax)]
    
    L_log_aRmax[1,1] = 1;
    L_log_aRmax[2,1] = rho_log_aRmax;
    L_log_aRmax[1,2] = 0;
    L_log_aRmax[2,2] = sqrt(1 - rho_log_aRmax^2);
    aRmax = append_col(log_a_z, log_Rmax_z);
    sigma_log_aRmax[1] = sigma_log_a;
    sigma_log_aRmax[2] = sigma_log_Rmax;
    aRmax = (diag_matrix(sigma_log_aRmax) * L_log_aRmax * aRmax')';
    a = exp(mu_log_a + col(aRmax,1));
    Rmax = exp(mu_log_Rmax + col(aRmax,2));
  }
  
  # AR(1) model for log(phi)
  phi[1] = log_phi_z[1]*sigma_log_phi/sqrt(1 - rho_log_phi^2); # initial anomaly
  for(i in 2:N_year)
    phi[i] = rho_log_phi*phi[i-1] + log_phi_z[i]*sigma_log_phi;
  phi = exp(phi);

  # Predict recruitment
  R_hat = rep_vector(0,N);
  for(i in 1:N_fit)
    R_hat[which_fit[i]] = A[which_fit[i]]*SR(a[pop[which_fit[i]]], Rmax[pop[which_fit[i]]], S[which_fit[i]], A[which_fit[i]]);
}

model {
  # Priors
  mu_log_a ~ normal(0,5);
  sigma_log_a ~ pexp(0,3,10);
  mu_log_Rmax ~ normal(0,10);
  sigma_log_Rmax ~ pexp(0,3,10);
  rho_log_phi ~ pexp(0,0.85,50);  # mildly regularize rho to ensure stationarity
  sigma_log_phi ~ pexp(0,3,10);
  sigma ~ pexp(0,2,10);
  
  # Hierarchical priors
  # [log(a), log(Rmax)] ~ MVN(0,D*R_log_aRmax*D), where D = diag_matrix(sigma_log_a, sigma_log_Rmax)
  log_a_z ~ normal(0,1);
  log_Rmax_z ~ normal(0,1);
  log_phi_z ~ normal(0,1);    # log(phi) ~ N(0, sigma_log_phi)
  
  # Likelihood
  R[which_fit] ~ lognormal(log(R_hat[which_fit]) + log(phi[year[which_fit]]), sigma);
}

generated quantities {
  vector[N] S_sim;    # simulated spawners
  vector[N] R_sim;    # simulated recruits

  S_sim = S;
  R_sim = R;

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

    if(R_NA[i] == 1)
      R_sim[i] = A[i] * SR(a[pop[i]], Rmax[pop[i]], S_sim[i], A[i]) * phi[year[i]] * lognormal_rng(0,sigma);
  }
}
