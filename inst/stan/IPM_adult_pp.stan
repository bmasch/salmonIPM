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
  
  # convert matrix to array of column vectors
  vector[] matrix_to_array(matrix m) {
    vector[2] arr[cols(m)];
    
    for(i in 1:cols(m))
      arr[i] = col(m,i);
    return(arr);
  }
  
  # Vectorized logical equality
  int[] veq(int[] x, int y) {
    int xeqy[size(x)];
    for(i in 1:size(x))
      xeqy[i] = x[i] == y;
    return(xeqy);
  }
  
  # Vectorized logical &&
  int[] vand(int[] cond1, int[] cond2) {
    int cond1_and_cond2[size(cond1)];
    for(i in 1:size(cond1))
      cond1_and_cond2[i] = cond1[i] && cond2[i];
    return(cond1_and_cond2);
  }
  
  # R-style conditional subsetting
  int[] rsub(int[] x, int[] cond) {
    int xsub[sum(cond)];
    int pos;
    pos = 1;
    for (i in 1:size(x))
      if (cond[i])
      {
        xsub[pos] = x[i];
        pos = pos + 1;
      }
    return(xsub);
  }
  
  # Equivalent of R: which(cond), where sum(cond) == 1
  int which(int[] cond) {
    int which_cond;
    for(i in 1:size(cond))
      if(cond[i])
        which_cond = i;
      return(which_cond);
  }
}

data {
  int<lower=1> N;                      # total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         # population identifier
  int<lower=1,upper=N> year[N];        # brood year identifier
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
  vector<lower=0>[N] A;                # habitat area associated with each spawner abundance obs
  vector<lower=0,upper=1>[N] F_rate;   # fishing mortality of wild adults
  int<lower=0,upper=N> N_B;            # number of years with B_take > 0
  int<lower=1,upper=N> which_B[max(N_B,1)]; # years with B_take > 0
  vector[max(N_B,1)] B_take_obs;       # observed broodstock take of wild adults
  int<lower=0> N_fwd;                  # total number of cases in forward simulations
  int<lower=1,upper=N> pop_fwd[max(N_fwd,1)]; # population identifier for forward simulations
  int<lower=1,upper=N+N_fwd> year_fwd[max(N_fwd,1)]; # brood year identifier for forward simulations
  vector<lower=0>[max(N_fwd,1)] A_fwd; # habitat area for each forward simulation
  vector<lower=0,upper=1>[max(N_fwd,1)] F_rate_fwd; # fishing mortality for forward simulations
  vector<lower=0,upper=1>[max(N_fwd,1)] B_rate_fwd; # broodstock take rate for forward simulations
  vector<lower=0,upper=1>[max(N_fwd,1)] p_HOS_fwd; # p_HOS for forward simulations
  int<lower=1> N_X;                    # number of productivity covariates
  matrix[max(max(year),max(year_fwd)),N_X] X; # brood-year productivity covariates (if none, use vector of zeros)
}

transformed data {
  int<lower=1,upper=N> N_pop;            # number of populations
  int<lower=1,upper=N> N_year;           # number of years, not including forward simulations
  int<lower=1,upper=N> N_year_all;       # total number of years, including forward simulations
  int<lower=2> ages[N_age];              # adult ages
  int<lower=0> n_HW_tot_obs[max(N_H,1)]; # total sample sizes for H/W frequencies
  int<lower=1> pop_year_indx[N];         # index of years within each pop, starting at 1
  int<lower=0,upper=N> fwd_init_indx[max(N_fwd,1),N_age]; # links "fitted" brood years to recruits in forward sims
  
  N_pop = max(pop);
  N_year = max(year);
  N_year_all = max(year_fwd);
  for(j in 1:N_age)
    ages[j] = max_age - N_age + j;
  for(i in 1:max(N_H,1)) n_HW_tot_obs[i] = n_H_obs[i] + n_W_obs[i];
  
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  
  fwd_init_indx = rep_array(0, max(N_fwd,1), N_age);
  if(N_fwd > 0)
  {
    for(i in 1:N_fwd)
    {
      for(j in 1:N_age)
      {
        if(year_fwd[i] - ages[j] < min(rsub(year_fwd, veq(pop_fwd, pop_fwd[i]))))
          fwd_init_indx[i,j] = which(vand(veq(pop, pop_fwd[i]), veq(year, year_fwd[i] - ages[j])));
      }
    }
  }
}

parameters {
  real mu_log_a;                        # hyper-mean log intrinsic productivity
  real<lower=0> sigma_log_a;            # hyper-SD log intrinsic productivity
  vector[N_pop] log_a_z;                # log intrinsic prod (Z-scores)
  real mu_log_Rmax;                     # hyper-mean log asymptotic recruitment
  real<lower=0> sigma_log_Rmax;         # hyper-SD log asymptotic recruitment
  vector[N_pop] log_Rmax_z;             # log asymptotic recruitment (Z-scores)
  real<lower=-1,upper=1> rho_log_aRmax; # correlation between log(a) and log(Rmax)
  vector[N_X] beta_log_phi;             # regression coefs for log productivity anomalies
  real<lower=-1,upper=1> rho_log_phi;   # AR(1) coef for log productivity anomalies
  real<lower=0> sigma_log_phi;          # hyper-SD of brood year log productivity anomalies
  vector[N_year] log_phi_z;             # log brood year productivity anomalies (Z-scores)
  real<lower=0> sigma_proc;             # unique process error SD
  simplex[N_age] mu_p;                  # mean of cohort age distributions
  vector<lower=0>[N_age-1] sigma_gamma; # among-pop SD of mean log-ratio age distributions
  cholesky_factor_corr[N_age-1] L_gamma; # Cholesky factor of among-pop correlation matrix of mean log-ratio age distns
  matrix[N_pop,N_age-1] gamma_z;        # population mean log-ratio age distributions (Z-scores)
  vector<lower=0>[N_age-1] sigma_alr_p; # SD of log-ratio cohort age distributions
  cholesky_factor_corr[N_age-1] L_alr_p; # Cholesky factor of correlation matrix of cohort log-ratio age distributions
  matrix[N,N_age-1] alr_p_z;            # log-ratio cohort age distributions (Z-scores)
  vector<lower=0>[max_age*N_pop] S_tot_init;  # true total spawner abundance in years 1-max_age
  simplex[N_age] q_init[max_age*N_pop]; # true wild spawner age distributions in years 1-max_age
  vector<lower=0,upper=1>[max(N_H,1)] p_HOS; # true p_HOS in years which_H
  vector[N] log_R_tot_z;                # log true recruit abundance (not density) by brood year (z-scores)
  vector<lower=0,upper=1>[max(N_B,1)] B_rate; # true broodstock take rate when B_take > 0
  real<lower=0> sigma_obs;              # observation error SD of total spawners
}

transformed parameters {
  vector<lower=0>[N_pop] a;           # intrinsic productivity 
  vector<lower=0>[N_pop] Rmax;        # asymptotic recruitment 
  vector<lower=0>[N_year] phi;        # brood year productivity anomalies
  vector<lower=0>[N] S_W_tot;         # true total wild spawner abundance
  vector[N] S_H_tot;                  # true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S_tot;           # true total spawner abundance
  row_vector[N_age-1] mu_alr_p;       # mean of log-ratio cohort age distributions
  matrix[N_pop,N_age-1] gamma;        # population mean log-ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p; # cohort age distributions
  matrix<lower=0,upper=1>[N,N_age] q; # true spawner age distributions
  vector[N] p_HOS_all;                # true p_HOS in all years (can == 0)
  vector<lower=0>[N] R_tot_hat;       # expected recruit abundance (not density) by brood year
  vector<lower=0>[N] R_tot;           # true recruit abundance (not density) by brood year
  vector<lower=0,upper=1>[N] B_rate_all; # true broodstock take rate in all years
  
  # Multivariate Matt trick for [log(a), log(b)]
  {
    matrix[2,2] L_log_aRmax;    # temp variable: Cholesky factor of correlation matrix of log(a) and log(Rmax)
    matrix[N_pop,2] aRmax;      # temp variable: matrix of a and Rmax
    vector[2] sigma_log_aRmax;  # temp variable: SD vector of [log(a), log(Rmax)]
    
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
  phi = exp(phi - mean(phi) + X[1:N_year,]*beta_log_phi);  # constrain log anomalies to sum to 0 (X should be centered)
  
  # Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  if(N_H > 0)
    p_HOS_all[which_H] = p_HOS;
  
  B_rate_all = rep_vector(0,N);
  if(N_B > 0)
    B_rate_all[which_B] = B_rate;
  
  # Multivariate Matt trick for age vectors (pop-specific mean and within-pop, time-varying)
  mu_alr_p = to_row_vector(log(mu_p[1:(N_age-1)]) - log(mu_p[N_age]));
  gamma = rep_matrix(mu_alr_p,N_pop) + (diag_matrix(sigma_gamma) * L_gamma * gamma_z')';
  p = append_col(gamma[pop,] + (diag_matrix(sigma_alr_p) * L_alr_p * alr_p_z')', rep_vector(0,N));
  
  # Calculate true total wild and hatchery spawners and spawner age distribution
  # and predict recruitment from brood year i
  for(i in 1:N)
  {
    row_vector[N_age] exp_p; # temp variable: exp(p[i,])
    row_vector[N_age] S_W;   # temp variable: true wild spawners by age
    
    # Inverse log-ratio transform of cohort age distn
    # (built-in softmax function doesn't accept row vectors)
    exp_p = exp(p[i,]);
    p[i,] = exp_p/sum(exp_p);
    
    if(pop_year_indx[i] <= max_age)
    {
      # Use initial values
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
    R_tot_hat[i] = A[i]*SR(a[pop[i]], Rmax[pop[i]], S_tot[i], A[i]);
    R_tot[i] = R_tot_hat[i]*phi[year[i]]*exp(sigma_proc*log_R_tot_z[i]);
  }
}

model {
  vector[max(N_B,1)] B_take; # true broodstock take when B_take_obs > 0
  
  # Priors
  mu_log_a ~ normal(2,5);
  sigma_log_a ~ pexp(0,3,10);
  mu_log_Rmax ~ normal(0,10);
  sigma_log_Rmax ~ pexp(0,3,10);
  beta_log_phi ~ normal(0,5);
  rho_log_phi ~ pexp(0,0.85,50);  # mildly regularize rho to ensure stationarity
  sigma_log_phi ~ pexp(0,2,10);
  sigma_proc ~ pexp(0,1,10);
  for(i in 1:(N_age-1))
  {
    sigma_gamma[i] ~ pexp(0,2,5);
    sigma_alr_p[i] ~ pexp(0,2,5); 
  }
  L_gamma ~ lkj_corr_cholesky(1);
  L_alr_p ~ lkj_corr_cholesky(1);
  sigma_obs ~ pexp(0,1,10);
  S_tot_init ~ lognormal(0,10);
  if(N_B > 0)
  {
    B_take = B_rate .* S_W_tot[which_B] .* (1 - q[which_B,1]) ./ (1 - B_rate);
    B_take_obs ~ lognormal(log(B_take), 0.1); # penalty to force pred and obs broodstock take to match 
  }
  
  # Hierarchical priors
  # [log(a), log(Rmax)] ~ MVN(0,D*R_log_aRmax*D), where D = diag_matrix(sigma_log_a, sigma_log_Rmax)
  log_a_z ~ normal(0,1);
  log_Rmax_z ~ normal(0,1);
  log_phi_z ~ normal(0,1);   # log(phi[i]) ~ N(rho_log_phi*log(phi[i-1]), sigma_log_phi)
  # pop mean age probs logistic MVN: gamma[i,] ~ MVN(0,D*R_gamma*D), where D = diag_matrix(sigma_gamma)
  to_vector(gamma_z) ~ normal(0,1);
  # age probs logistic MVN: alr_p[i,] ~ MVN(gamma[pop[i],], D*R_alr_p*D), where D = diag_matrix(sigma_alr_p)
  to_vector(alr_p_z) ~ normal(0,1);
  
  # Process model
  log_R_tot_z ~ normal(0,1); # total recruits: R_tot ~ lognormal(log(R_tot_hat), sigma_proc)
  
  # Observation model
  S_tot_obs[which_S_obs] ~ lognormal(log(S_tot[which_S_obs]), sigma_obs);   # observed total spawners
  if(N_H > 0) n_H_obs ~ binomial(n_HW_tot_obs, p_HOS); # observed counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));                  # obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  corr_matrix[N_age-1] R_gamma;               # among-pop correlation matrix of mean log-ratio age distns 
  corr_matrix[N_age-1] R_alr_p;               # correlation matrix of within-pop cohort log-ratio age distns 
  # vector<lower=0>[N_year_all] phi_all;        # brood year productivity anomalies in all years
  # vector<lower=0>[N_fwd] S_W_tot_fwd;         # true total wild spawner abundance in forward simulations
  # vector[N_fwd] S_H_tot_fwd;                  # true total hatchery spawner abundance in forward simulations
  # vector<lower=0>[N_fwd] S_tot_fwd;           # true total spawner abundance in forward simulations
  # matrix<lower=0,upper=1>[N_fwd,N_age] p_fwd; # cohort age distributions in forward simulations
  # matrix<lower=0,upper=1>[N_fwd,N_age] q_fwd; # spawner age distributions in forward simulations
  # vector<lower=0>[N_fwd] R_tot_hat_fwd;       # expected recruit abundance by brood year in forward simulations
  # vector<lower=0>[N_fwd] R_tot_fwd;           # true recruit abundance by brood year in forward simulations
  # vector[N] ll_S_tot_obs;                     # pointwise log-likelihood of total spawners
  # vector[N_H] ll_n_H_obs;                     # pointwise log-likelihood of hatchery vs. wild frequencies
  # vector[N] ll_n_age_obs;                     # pointwise log-likelihood of wild age frequencies
  
  R_gamma = multiply_lower_tri_self_transpose(L_gamma);
  R_alr_p = multiply_lower_tri_self_transpose(L_alr_p);
  
  # # AR(1) model for log(phi)
  # phi_all[1:N_year] = phi; # copy from phi
  # for(i in (N_year + 1):N_year_all)
  #   phi_all[i] = lognormal_rng(rho_log_phi*log(phi_all[i-1]) + X[i,]*beta_log_phi, sigma_log_phi);
  # 
  # # Calculate true total wild and hatchery spawners and spawner age distribution
  # # and simulate recruitment from brood year i
  # # (Note that if N_fwd == 0, this block will not execute)
  # for(i in 1:N_fwd)
  # {
  #   vector[N_age-1] alr_p_fwd;   # temp variable: alr(p_fwd[i,])'
  #   row_vector[N_age] S_W_fwd;   # temp variable: true wild spawners by age
  #   
  #   # Inverse log-ratio transform of cohort age distn
  #   alr_p_fwd = multi_normal_cholesky_rng(to_vector(gamma[pop_fwd[i],]), L_alr_p);
  #   p_fwd[i,] = to_row_vector(softmax(append_row(alr_p_fwd,0)));
  #   
  #   for(j in 1:N_age)
  #   {
  #     if(fwd_init_indx[i,j] != 0)
  #     {
  #       # Use estimated values from previous cohorts
  #       S_W_fwd[j] = R_tot[fwd_init_indx[i,j]]*p[fwd_init_indx[i,j],j];
  #     }
  #     else
  #     {
  #       S_W_fwd[j] = R_tot_fwd[i-ages[j]]*p_fwd[i-ages[j],j];
  #     }
  #   }
  #   
  #   for(j in 2:N_age)  # catch and broodstock removal (assumes no take of age 1)
  #     S_W_fwd[j] = S_W_fwd[j]*(1 - F_rate_fwd[i])*(1 - B_rate_fwd[i]);
  #   
  #   S_W_tot_fwd[i] = sum(S_W_fwd);
  #   S_H_tot_fwd[i] = S_W_tot_fwd[i]*p_HOS_fwd[i]/(1 - p_HOS_fwd[i]);
  #   q_fwd[i,] = S_W_fwd/S_W_tot_fwd[i];
  #   S_tot_fwd[i] = S_W_tot_fwd[i] + S_H_tot_fwd[i];
  #   R_tot_hat_fwd[i] = A_fwd[i]*SR(a[pop_fwd[i]], Rmax[pop_fwd[i]], S_tot_fwd[i], A_fwd[i]);
  #   R_tot_fwd[i] = lognormal_rng(log(R_tot_hat_fwd[i]) + log(phi_all[year_fwd[i]]), sigma_proc);
  # }
  
  # ll_S_tot_obs = rep_vector(0,N);
  # for(i in 1:N_S_obs)
  #   ll_S_tot_obs[which_S_obs[i]] = lognormal_lpdf(S_tot_obs[which_S_obs[i]],
  #                                                 log(S_tot[which_S_obs[i]]), sigma_obs);
  # 
  # if(N_H > 0)
  # {
  #   for(i in 1:N_H)
  #     ll_n_H_obs[i] = binomial_lpmf(n_H_obs[i], n_HW_tot_obs[i], p_HOS[i]);
  # }
  # 
  # ll_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
}
