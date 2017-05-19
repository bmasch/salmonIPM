functions {
  # spawner-recruit functions
  real SR(real a, real b, real S, real A) {
    real R;
    R = a*S/(A + b*S);
    return(R);
  }
  
  # Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(real y, real mu, real sigma, real shape) {
    return(-(fabs(y - mu)/sigma)^shape);
  }
}

data {
  int<lower=1> N;                      # total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         # population identifier
  int<lower=1,upper=N> year[N];        # brood year identifier
  int<lower=1> N_X;                    # number of productivity covariates
  matrix[max(year),N_X] X;             # brood-year productivity covariates (if none, use vector of zeros)
  int<lower=0,upper=max(pop)> N_pop_H; # number of populations with hatchery input
  int<lower=1,upper=max(pop)> which_pop_H[max(N_pop_H,1)]; # populations with hatchery input
  int<lower=1,upper=N> N_S_obs;        # number of cases with non-missing spawner abundance obs 
  int<lower=1,upper=N> which_S_obs[N_S_obs]; # cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_tot_obs;        # observed annual total spawner abundance (not density)
  int<lower=2> N_age;                  # number of adult age classes
  int<lower=2> max_age;                # maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  # observed wild spawner age frequencies (all zero row = NA)  
  int<lower=0,upper=N> N_H;            # number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[max(N_H,1)]; # years with p_HOS > 0
  int<lower=0> n_W_obs[max(N_H,1)];    # count of wild spawners in samples (assumes no NAs)
  int<lower=0> n_H_obs[max(N_H,1)];    # count of hatchery spawners in samples (assumes no NAs)
  vector[N] A;                         # habitat area associated with each spawner abundance obs
  vector[N] F_rate;                    # fishing mortality rate of wild adults (no fishing on jacks)
  int<lower=0,upper=N> N_B;            # number of years with B_take > 0
  int<lower=1,upper=N> which_B[max(N_B,1)]; # years with B_take > 0
  vector[max(N_B,1)] B_take_obs;       # observed broodstock take of wild adults
}

transformed data {
  int<lower=1,upper=N> N_pop;            # number of populations
  int<lower=1,upper=N> N_year;           # number of years
  int<lower=2> ages[N_age];              # adult ages
  int<lower=1> pop_year_indx[N];         # index of years within each pop, starting at 1
  int<lower=0> n_HW_tot_obs[max(N_H,1)]; # total sample sizes for H/W frequencies
  
  N_pop = max(pop);
  N_year = max(year);
  for(j in 1:N_age)
    ages[j] = max_age - N_age + j;
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  for(i in 1:max(N_H,1)) n_HW_tot_obs[i] = n_H_obs[i] + n_W_obs[i];
}

parameters {
  vector<lower=0>[N_pop] a;             # intrinsic productivity
  vector<lower=0>[N_pop] b;             # density dependence
  matrix[N_pop,N_X] beta_proc;          # regression coefs for log productivity anomalies
  vector<lower=-1,upper=1>[N_pop] rho_proc; # AR(1) coefs for log productivity anomalies
  vector<lower=0>[N_pop] sigma_proc;    # process error SDs
  simplex[N_age] exp_gamma[N_pop];      # population mean age distributions
  matrix<lower=0>[N_pop,N_age-1] sigma_alr_p; # log-ratio cohort age distribution SDs
  cholesky_factor_corr[N_age-1] L_alr_p[N_pop]; # Cholesky factors of correlation matrices of cohort log-ratio age distributions
  matrix[N,N_age-1] alr_p_z;            # log-ratio cohort age distributions (Z-scores)
  vector<lower=0>[max_age*N_pop] S_tot_init;  # true total spawner abundance in years 1-max_age
  simplex[N_age] q_init[max_age*N_pop]; # true wild spawner age distributions in years 1-max_age
  vector<lower=0,upper=1>[max(N_H,1)] p_HOS; # true p_HOS in years which_H
  vector[N] log_R_tot_z;                # log true recruit abundance (not density) by brood year (z-scores)
  vector<lower=0,upper=1>[max(N_B,1)] B_rate; # true broodstock take rate when B_take > 0
  vector<lower=0>[N_pop] sigma_obs;     # observation error SDs
}

transformed parameters {
  vector<lower=0>[N] S_W_tot;         # true total wild spawner abundance
  vector[N] S_H_tot;                  # true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S_tot;           # true total spawner abundance
  matrix[N_pop,N_age-1] gamma;        # population mean log ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p; # cohort age distributions
  matrix<lower=0,upper=1>[N,N_age] q; # true spawner age distributions
  vector[N] p_HOS_all;                # true p_HOS in all years (can == 0)
  vector<lower=0>[N] R_tot_hat;       # expected recruit abundance (not density) by brood year
  vector[N] log_R_tot_proc;           # process error in recruit abundance by brood year 
  vector<lower=0>[N] R_tot;           # true recruit abundance (not density) by brood year
  vector<lower=0,upper=1>[N] B_rate_all; # true broodstock take rate in all years
  
  # Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  if(N_H > 0)
    p_HOS_all[which_H] = p_HOS;
  
  B_rate_all = rep_vector(0,N);
  if(N_B > 0)
    B_rate_all[which_B] = B_rate;
  
  # Multivariate Matt trick for within-pop, time-varying age vectors
  for(j in 1:N_pop)
    gamma[j,] = to_row_vector(log(exp_gamma[j,1:(N_age-1)]) - log(exp_gamma[j,N_age]));

  # Calculate true total wild and hatchery spawners and spawner age distribution
  # and predict recruitment from brood year t
  for(i in 1:N)
  {
    row_vector[N_age] alr_p; # temp variable: alr(p[i,])
    row_vector[N_age] S_W;   # temp variable: true wild spawners by age
    
    # inverse log-ratio transform of cohort age distn
    # (built-in softmax function doesn't accept row vectors)
    alr_p = rep_row_vector(0,N_age);
    alr_p[1:(N_age-1)] = gamma[pop[i],] + to_row_vector(diag_matrix(to_vector(sigma_alr_p[pop[i],])) * L_alr_p[pop[i]] * to_vector(alr_p_z[i,]));
    alr_p = exp(alr_p);
    p[i,] = alr_p/sum(alr_p);
    
    if(pop_year_indx[i] <= max_age)
    {
      # use initial values
      S_W_tot[i] = S_tot_init[(pop[i]-1)*max_age+pop_year_indx[i]]*(1 - p_HOS_all[i]);        
      S_H_tot[i] = S_tot_init[(pop[i]-1)*max_age+pop_year_indx[i]]*p_HOS_all[i];
      q[i,1:N_age] = to_row_vector(q_init[(pop[i]-1)*max_age+pop_year_indx[i],1:N_age]);
      S_W = S_W_tot[i]*q[i,];
    }
    else
    {
      for(j in 1:N_age)
        S_W[j] = R_tot[i-ages[j]]*p[i-ages[j],j];
      for(j in 2:N_age)  # catch and broodstock removal (assumes no take of age 1)
        S_W[j] = S_W[j]*(1 - F_rate[i])*(1 - B_rate_all[i]);
      S_W_tot[i] = sum(S_W);
      S_H_tot[i] = S_W_tot[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
      q[i,] = S_W/S_W_tot[i];
    }
    
    S_tot[i] = S_W_tot[i] + S_H_tot[i];
    R_tot_hat[i] = A[i]*SR(a[pop[i]], b[pop[i]], S_tot[i], A[i]);
    if(pop_year_indx[i] == 1) # initial process error
      log_R_tot_proc[i] = log_R_tot_z[i]*sigma_proc[pop[i]];
    else
      log_R_tot_proc[i] = rho_proc[pop[i]]*log_R_tot_proc[i-1] + log_R_tot_z[i]*sigma_proc[pop[i]];
    log_R_tot_proc[i] =  dot_product(X[year[i],], beta_proc[pop[i],]) + log_R_tot_proc[i];
    R_tot[i] = R_tot_hat[i]*exp(log_R_tot_proc[i]);
  }
}

model {
  vector[max(N_B,1)] B_take; # true broodstock take when B_take_obs > 0
  
  # Priors
  a ~ normal(0,500);
  b ~ normal(0,1000);
  to_vector(beta_proc) ~ normal(0,5);
  for(j in 1:N_pop)
  {
    rho_proc[j] ~ pexp(0,0.85,20); # mildly regularize rho to ensure stationarity
    L_alr_p[j] ~ lkj_corr_cholesky(3);
    sigma_proc[j] ~ pexp(0,1,10);
    sigma_obs[j] ~ pexp(1,0.85,20);  # avoid sigma_obs < 0.1, which leads to divergences
  }
  to_vector(sigma_alr_p) ~ normal(0,5);
  S_tot_init ~ lognormal(0,5);
  if(N_B > 0)
  {
    B_take = B_rate .* S_W_tot[which_B] .* (1 - q[which_B,1]) ./ (1 - B_rate);
    B_take_obs ~ lognormal(log(B_take), 0.1); # penalty to force pred and obs broodstock take to match 
  }
  
  # Hierarchical priors
  # age probs logistic MVN: alr_p[i,] ~ MVN(gamma[pop[i],], D*R_alr_p*D), where D = diag_matrix(sigma_alr_p)
  to_vector(alr_p_z) ~ normal(0,1);
  
  # Process model
  log_R_tot_z ~ normal(0,1); # total recruits: R_tot ~ lognormal(log(R_tot_hat), sigma_proc)
  
  # Observation model
  S_tot_obs[which_S_obs] ~ lognormal(log(S_tot[which_S_obs]), sigma_obs[pop[which_S_obs]]);   # observed total spawners
  if(N_H > 0) n_H_obs ~ binomial(n_HW_tot_obs, p_HOS); # observed counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));                 # obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  corr_matrix[N_age-1] R_alr_p[N_pop]; # correlation matrices of within-pop cohort log-ratio age distns
  # vector[N] ll_S_tot_obs;  # pointwise log-likelihood of total spawners
  # vector[N_H] ll_n_H_obs;  # pointwise log-likelihood of hatchery vs. wild frequencies
  # vector[N] ll_n_age_obs;  # pointwise log-likelihood of wild age frequencies

  for(j in 1:N_pop)
    R_alr_p[j] = multiply_lower_tri_self_transpose(L_alr_p[j]);
  # ll_S_tot_obs = rep_vector(0,N);
  # for(i in 1:N_S_obs)
  #   ll_S_tot_obs[which_S_obs[i]] = lognormal_lpdf(S_tot_obs[which_S_obs[i]],
  #                                                 log(S_tot[which_S_obs[i]]),
  #                                                 sigma_obs[pop[which_S_obs[i]]]);
  # 
  # if(N_H > 0)
  # {
  #   for(i in 1:N_H)
  #     ll_n_H_obs[i] = binomial_lpmf(n_H_obs[i], n_HW_tot_obs[i], p_HOS[i]);
  # }
  # 
  # ll_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
}
