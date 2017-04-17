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
  int<lower=1,upper=N> N_pop;          # number of populations
  int<lower=1,upper=N> N_year;         # number of years
  int<lower=1> pop_indx[max(pop),2];   # start and end position of each pop
  
  N_pop = max(pop);
  N_year = max(year);
  pop_indx[pop[1],1] = 1;
  for(i in 2:N)
  {
    if(pop[i] != pop[i-1])
    {
      pop_indx[pop[i-1],2] = i - 1;
      pop_indx[pop[i],1] = 1;
    }
    pop_indx[pop[N],2] = N;
  }
}

parameters {
  vector<lower=0>[N_pop] a;             # intrinsic productivity
  vector<lower=0>[N_pop] b;             # density dependence
  vector<lower=-1,upper=1>[N_pop] rho;  # AR(1) coefs of residuals
  vector<lower=0>[N_pop] sigma;         # residual error SD
}

transformed parameters {
  vector<lower=0>[N] R_hat;             # expected recruit abundance (not density) by brood year
  
  # Predict recruitment
  R_hat = A .* SR(a[pop], b[pop], S, A);
}

model {
  # Priors
  a ~ lognormal(0,5);
  b ~ lognormal(0,5);
  for(i in 1:N_pop)
    sigma[i] ~ pexp(0,2,10);

  # Likelihood
  for(j in 1:N_pop)
  {
    int N_j;
    vector[pop_indx[j,2]-pop_indx[j,1]+1] err_j;
    
    err_j = log(R[pop_indx[j,1]:pop_indx[j,2]]) - log(R_hat[pop_indx[j,1]:pop_indx[j,2]]);
    N_j = num_elements(err_j);
    tail(err_j, N_j - 1) ~ normal(rho[j]*head(err_j, N_j - 1), sigma[j]);
  }
}
