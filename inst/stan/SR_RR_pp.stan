functions {
  # spawner-recruit functions
  vector SR(vector a, vector b, vector S, vector A) {
    vector[num_elements(S)] R;
    R = a .* S ./ (A + b .* S);
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
  vector<lower=0>[N] S;         # observed annual total spawner abundance (not density)
  vector<lower=0>[N] R;         # total natural recruit abundance (not density), including harvest and broodstock removals
  vector[N] A;                  # habitat area associated with each spawner abundance obs
}

transformed data {
  int<lower=1,upper=N> N_pop;            # number of populations
  int<lower=1,upper=N> N_year;           # number of years
  
  N_pop = max(pop);
  N_year = max(year);
}

parameters {
  real mu_log_a;                        # hyper-mean log intrinsic productivity of wild spawners
  real<lower=0> sigma_log_a;            # hyper-SD log intrinsic productivity
  vector[N_pop] log_a_z;                # log intrinsic prod of wild spawners (Z-scores)
  real mu_log_b;                        # hyper-mean log density dependence
  real<lower=0> sigma_log_b;            # hyper-SD log density dependence
  vector[N_pop] log_b_z;                # log density dependence (Z-scores)
  real<lower=-1,upper=1> rho_log_phi;   # AR(1) coef for brood year log productivity anomalies
  real<lower=0> sigma_log_phi;          # hyper-SD of brood year log productivity anomalies
  vector[max(year)] log_phi_z;          # log brood year productivity anomalies (Z-scores)
  real<lower=0> sigma;                  # residual error SD
}

transformed parameters {
  vector<lower=0>[N_pop] a;             # intrinsic productivity of wild spawners
  vector<lower=0>[N_pop] b;             # density dependence 
   vector<lower=0>[N_year] phi;         # brood year productivity anomalies
  vector<lower=0>[N] R_hat;             # expected recruit abundance (not density) by brood year
  
  phi[1] = log_phi_z[1]*sigma_log_phi; # initial anomaly
  for(i in 2:N_year)
    phi[i] = rho_log_phi*phi[i-1] + log_phi_z[i]*sigma_log_phi;
  phi = exp(phi);
  a = exp(mu_log_a + sigma_log_a*log_a_z);
  b = exp(mu_log_b + sigma_log_b*log_b_z);

  # Predict recruitment
  R_hat = A .* SR(a[pop], b[pop], S, A) .* phi[year];
}

model {
  # Priors
  mu_log_a ~ normal(0,5);
  sigma_log_a ~ pexp(0,3,10);
  mu_log_b ~ normal(0,5);
  sigma_log_b ~ pexp(0,3,10);
  # rho_log_phi ~ pexp(0,0.8,10);  # mildly regularize rho to ensure stationarity
  sigma_log_phi ~ pexp(0,1,10);
  sigma ~ pexp(0,1,10);

  # Hierarchical priors
  log_a_z ~ normal(0,1);      # log(a) ~ N(mu_log_a, sigma_log_a)
  log_b_z ~ normal(0,1);      # log(b) ~ N(mu_log_b, sigma_log_b)
  log_phi_z ~ normal(0,1);    # log(phi) ~ N(0, sigma_log_phi)
  
  # Likelihood
  R ~ lognormal(log(R_hat), sigma);
}
