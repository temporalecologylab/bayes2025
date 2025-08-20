// See Mike's repo on HMM: https://github.com/betanalpha/quarto_modeling_techniques/tree/main/hidden_markov_models

functions {
  // Mean-dispersion parameterization of negative-binomial family
  // 0 < mu < +inf
  // 0 < psi    < +inf
  real neg_binomial_alt_lpmf(int n, real mu, real psi) {
    return neg_binomial_2_lpmf(n | mu, 1 / psi);
  }

  int neg_binomial_2_alt_log_safe_rng(real mu, real psi) {
    real psihere = 1/psi;
    real gamma_rate = gamma_rng(psihere, psihere / exp(mu));
    if (gamma_rate > exp(20.7)) gamma_rate = exp(20.7); // max value before overflow
    return poisson_rng(gamma_rate);
  }
  
}

data {
  int<lower=1> N;           // Total number of observations
  int<lower=1> N_trees;     // Number of trees
  int<lower=1> N_max_years; // total number of years 

  // Number of observations per tree
  array[N_trees] int<lower=1, upper=N> N_years;

  // Ragged array indexing
  array[N_trees] int<lower=1, upper=N> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N> tree_end_idxs;

  // Number of seeds in each tree each year
  array[N] int<lower=0> seed_counts;
  
  // Tree observation years
  array[N] int<lower=1, upper=N_max_years> years;
  
  // Previous summer temperature (in degC)
  array[N_max_years*N_trees] real prevsummer_temps;
  
  // GDD to last spring frost (in x10 degC)
  array[N_max_years*N_trees] real gdd_lastfrost;
  
  // Spring temperature (in degC)
  array[N_max_years*N_trees] real spring_temps;
  
  // below, for predictions!
  int<lower=1> Nnew; 
  int<lower=1> N_newtrees; 
  int<lower=1> N_max_newyears; 
  array[N_newtrees] int<lower=1, upper=Nnew> N_newyears;
  array[N_newtrees] int<lower=1, upper=Nnew> newtree_start_idxs;
  array[N_newtrees] int<lower=1, upper=Nnew> newtree_end_idxs;
  array[Nnew] int<lower=1, upper=N_max_newyears> newyears;
  array[N_newtrees] int<lower=1> newtree_stand_idxs;
  array[Nnew] real newprevsummer_temps;
  array[Nnew] real newgdd_lastfrost;
  array[Nnew] real newspring_temps;
}

transformed data{
   int N_newstands = Nnew/N_newtrees;
}

parameters {
  real<lower=0> lambda1; // Non-masting intensity
  real<lower=0> psi1; // Non-masting dispersion
  real<lower=0, upper=1> theta1; // probability of drawing a zero (zero-inflation)
  
  real<lower=lambda1> lambda20; // Masting intensity
  real<upper=0> beta_lambda2_frost; // Effect of frost risk on masting intensity
  real<lower=0> beta_lambda2_spring; // Effect of spring temp. on masting intensity
  real<lower=0> psi2; // Masting dispersion

  real<lower=0, upper=1> rho0; // Initial masting probability
  real<lower=0, upper=1> tau_nm_m0; // Non-masting to masting probability
  real<lower=0> beta_nm_m; // Effect of previous summer temp. on transition prob. to masting 
  real<lower=0, upper=1> tau_m_m0; // Masting to masting probability
  
}

transformed parameters {
  
  simplex[2] rho = [1 - rho0, rho0]';
  
  // define baseline conditions
  real summertemp0 = 15;
  real frostgdd0 = 15;
  real springtemp0 = 7;
  
  real nugget = 1e-20;
  
}

model {

  lambda1 ~ normal(0, 20 / 2.57); 
  psi1 ~ normal(0, 5 / 2.57); 
  lambda20 ~ normal(0, 500 / 2.57); 
  beta_lambda2_frost ~ normal(0, log(1.2) / 2.57); # 20% relative change with one unit of frost risk
  beta_lambda2_spring ~ normal(0, log(1.2) / 2.57); # 20% relative change with one unit of frost risk
  psi2 ~ normal(0, 5 / 2.57); 
  beta_nm_m ~ normal(0, 0.6 / 2.57); # ~20% relative change with one unit of frost risk, at tau0 = 0.5
  
  theta1 ~ beta(2, 2); 
  tau_nm_m0 ~ beta(1, 1.5); 
  tau_m_m0 ~ beta(1, 1.5); 
  // Implicit uniform prior model over rho0

  for (t in 1:N_trees) {
    
    int N_years_tree = N_years[t];
    array[N_years_tree] int seed_counts_tree = seed_counts[tree_start_idxs[t]:tree_end_idxs[t]];
    array[N_years_tree] int years_tree = years[tree_start_idxs[t]:tree_end_idxs[t]];
    array[N_max_years] real prevsummer_temps_tree = prevsummer_temps[(1+N_max_years*(t-1)):(N_max_years*t)];
    array[N_max_years] real gdd_lastfrost_tree = gdd_lastfrost[(1+N_max_years*(t-1)):(N_max_years*t)];
    array[N_max_years] real spring_temps_tree = spring_temps[(1+N_max_years*(t-1)):(N_max_years*t)];
    
    // Observational model
    {
      // Forward algorithm
      real norm;
      real log_norm = 0;

      vector[2] alpha = rho;

      norm = max(alpha);
      log_norm += log(norm);
      alpha /= norm;

      for (n in 1:N_years_tree) {
        
        int y = seed_counts_tree[n];
        
        vector[2] log_omega;
        // Non-masting
        if (y == 0){
          log_omega[1] = log_sum_exp(bernoulli_lpmf(1 | theta1), bernoulli_lpmf(0 | theta1) 
          + neg_binomial_alt_lpmf(y | lambda1, psi1));
        }else{
          log_omega[1] = bernoulli_lpmf(0 | theta1) + neg_binomial_alt_lpmf(y | lambda1, psi1);
        }
        // Masting
        real lambda2 = exp(log(lambda20) 
        + beta_lambda2_frost * (gdd_lastfrost_tree[years_tree[n]]-frostgdd0) 
        + beta_lambda2_spring * (spring_temps_tree[years_tree[n]]-springtemp0));
        log_omega[2] = neg_binomial_alt_lpmf(y | lambda2, psi2); 
        
        vector[2] omega;
        omega[1] = exp(log_omega[1]);
        omega[2] = exp(log_omega[2]);
        
        if(n == 1){
          if(years_tree[n] > 1){
            for(i in 2:years_tree[n]){
              // Construct transition matrix
              real tau_nm_m = inv_logit(logit(tau_nm_m0) + beta_nm_m * (prevsummer_temps_tree[i]-summertemp0));
              real tau_m_m = tau_m_m0;
              matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                                   [1 - tau_m_m, tau_m_m] ];
              alpha = Gamma' * alpha;
            }
          }
          alpha = omega .* alpha;
        }else{
          int delta = years_tree[n] - years_tree[n - 1];
          for (d in 1:delta){
            // Construct transition matrix
            real tau_nm_m = inv_logit(logit(tau_nm_m0) + beta_nm_m * (prevsummer_temps_tree[years_tree[n-1]+d]-summertemp0));
            real tau_m_m = tau_m_m0;
            matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                                 [1 - tau_m_m, tau_m_m] ];
            alpha = Gamma' * alpha;
          }
          alpha = omega .* alpha;
        }

        norm = max(alpha);
        log_norm += log(norm);
        alpha /= norm;
      }

      // Marginal observational model
      target += sum(alpha) + log_norm;
    }
    
  }
  
}

generated quantities {
  array[N_max_years*N_trees] int states_pred;       // Posterior latent states
  array[N_max_years*N_trees] real seed_counts_pred; // Posterior predictive observations
  array[N_max_years*N_trees] real tau_nm_m_pred; // Posterior predictive observations

  for (t in 1:N_trees) {

    int N_years_tree = N_years[t];
    array[N_years_tree] int seed_counts_tree = seed_counts[tree_start_idxs[t]:tree_end_idxs[t]];
    array[N_years_tree] int years_tree = years[tree_start_idxs[t]:tree_end_idxs[t]];
    int global_start_idx = (1+N_max_years*(t-1));
    int global_end_idx = (N_max_years*t);
    array[N_max_years] real prevsummer_temps_tree = prevsummer_temps[global_start_idx:global_end_idx];
    array[N_max_years] real gdd_lastfrost_tree = gdd_lastfrost[global_start_idx:global_end_idx];
    array[N_max_years] real spring_temps_tree = spring_temps[global_start_idx:global_end_idx];

    matrix[2, N_max_years] omega = rep_matrix(1, 2, N_max_years);

    
    // Construct transition matrix
    // matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
    //                      [tau_m_nm, 1 - tau_m_nm] ];

    // Forward algorithm

    // Build omega with the observed latent states
    for(n in 1:N_years_tree){
      vector[2] log_omega;
      int idx = years_tree[n];
      int y = seed_counts_tree[n];

      // Non-masting
      if (y == 0){
        log_omega[1] = log_sum_exp(bernoulli_lpmf(1 | theta1), bernoulli_lpmf(0 | theta1)
        + neg_binomial_2_log_lpmf(y | log(lambda1), 1/psi1));
      }else{
        log_omega[1] = bernoulli_lpmf(0 | theta1) + neg_binomial_2_log_lpmf(y | log(lambda1), 1/psi1);  // neg_binomial_alt_lpmf(y | lambda1, psi1)
      }
      // Masting
      real loglambda2 = log(lambda20) 
      + beta_lambda2_frost * (gdd_lastfrost_tree[years_tree[n]]-frostgdd0) 
      + beta_lambda2_spring * (spring_temps_tree[years_tree[n]]-springtemp0);
      
      log_omega[2] = neg_binomial_2_log_lpmf(y | loglambda2, 1/psi2); // neg_binomial_alt_lpmf(y | lambda2, psi2)

      omega[1,idx] = exp(log_omega[1])+nugget; // modified this, adding a small nugget
      omega[2,idx] = exp(log_omega[2])+nugget; // modified this, adding a small nugget
    }

    // Loop over all years
    array[N_max_years] vector[2] alpha;
    for (n in 1:N_max_years) {
      int global_n = global_start_idx + n - 1;

      real tau_nm_m = inv_logit(logit(tau_nm_m0) + beta_nm_m * (prevsummer_temps_tree[n]-summertemp0));
      real tau_m_m = tau_m_m0;
      matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                           [1 - tau_m_m, tau_m_m] ];
      tau_nm_m_pred[global_n] = tau_nm_m;

      if(n == 1){
        alpha[n] = omega[,n] .* rho;
      }
      else{
        alpha[n] = omega[,n] .* (Gamma' * alpha[n - 1]);
      }
      alpha[n] /= max(alpha[n]);
    }

    // Sample final latent state (p. 18 in Mike's HMM chapter)
    vector[2] r = alpha[N_max_years];
    vector[2] lambda = r / sum(r);
    // print(r);
    // print(sum(r));
    // print(lambda);
    // print("----------------");

    states_pred[global_end_idx] = categorical_rng(lambda);
    if(states_pred[global_end_idx] == 1) {
      if(bernoulli_rng(theta1)){
        seed_counts_pred[global_end_idx] = 0;
      }else{
        seed_counts_pred[global_end_idx] = neg_binomial_2_alt_log_safe_rng(log(lambda1), psi1);  //neg_binomial_alt_rng(lambda1, psi1)
      }
    }
    else if(states_pred[global_end_idx] == 2) {
      real loglambda2 = log(lambda20) 
      + beta_lambda2_frost * (gdd_lastfrost_tree[N_max_years]-frostgdd0) 
      + beta_lambda2_spring * (spring_temps_tree[N_max_years]-springtemp0);
      
      seed_counts_pred[global_end_idx] = neg_binomial_2_alt_log_safe_rng(loglambda2, psi2); // neg_binomial_alt_rng(lambda2, psi2) 
    }

    // Sample latent states while running the backward algorithm
    vector[2] beta = ones_vector(2);

    for (rn in 2:N_max_years) {
      int n = N_max_years + 1 - rn;
      int global_n = global_end_idx + 1 - rn;
      int states_pred_prev = states_pred[n + 1];

      vector[2] omega_prev = omega[, n + 1];

      real tau_nm_m = inv_logit(logit(tau_nm_m0) + beta_nm_m * (prevsummer_temps_tree[n+1]-summertemp0)); // was prevsummer_temps_tree[n] before
      real tau_m_m = tau_m_m0;
      matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                           [1 - tau_m_m, tau_m_m] ];

      r =   beta[states_pred_prev] * omega_prev[states_pred_prev]
          * ( alpha[n] .* Gamma[,states_pred_prev] );
      lambda = r / sum(r);
      
      // print(n);
      // print(beta);
      // print(omega_prev);
      // print(alpha[n]);
      // print(Gamma[,states_pred_prev]);
      // print("---");
      // print(r);
      // print(sum(r));
      // print(lambda);
      // print("----------------");

      states_pred[global_n] = categorical_rng(lambda);

      if(states_pred[global_n] == 1) {
        if(bernoulli_rng(theta1)){
          seed_counts_pred[global_n] = 0;
        }else{
          seed_counts_pred[global_n] = neg_binomial_2_alt_log_safe_rng(log(lambda1), psi1);  //neg_binomial_alt_rng(lambda1, psi1)
        }
      }
      else if(states_pred[global_n] == 2) {
        real loglambda2 = log(lambda20) 
        + beta_lambda2_frost * (gdd_lastfrost_tree[n]-frostgdd0) 
        + beta_lambda2_spring * (spring_temps_tree[n]-springtemp0);
        
        seed_counts_pred[global_n] = neg_binomial_2_alt_log_safe_rng(loglambda2, psi2); // neg_binomial_alt_rng(lambda2, psi2)
      }

      tau_nm_m = inv_logit(logit(tau_nm_m0) + beta_nm_m * (prevsummer_temps_tree[n+1]-summertemp0));
      tau_m_m = tau_m_m0;
      Gamma = [ [1 - tau_nm_m, tau_nm_m],
              [1 - tau_m_m, tau_m_m] ];
      beta = Gamma * (omega_prev .* beta); // p. 16 in Mike's HMM chapter
      beta /= max(beta);

      
    }
  }
  
  // Predict for new (fake) trees!
  array[Nnew] int states_new;
  array[Nnew] int seed_counts_new;
  
  array[N_newstands, N_max_newyears] int states_new_plot = rep_array(0, N_newstands, N_max_newyears);
  array[N_newstands, N_max_newyears] int seed_counts_new_plot = rep_array(0, N_newstands, N_max_newyears);
  
  // int global_start_idx = 1;
  // int global_end_idx = N_max_years;
  // array[N_max_years] real prevsummer_temps_tree = prevsummer_temps[global_start_idx:global_end_idx];
  // array[N_max_years] real gdd_lastfrost_tree = gdd_lastfrost[global_start_idx:global_end_idx];
  // array[N_max_years] real spring_temps_tree = spring_temps[global_start_idx:global_end_idx];
  
  
  for(newt in 1:N_newtrees){
    
    matrix[2, N_max_newyears] omega = rep_matrix(1, 2, N_max_newyears);
    
    array[N_max_newyears] real prevsummer_temps_tree = newprevsummer_temps[newtree_start_idxs[newt]:newtree_end_idxs[newt]];
    array[N_max_newyears] real gdd_lastfrost_tree = newgdd_lastfrost[newtree_start_idxs[newt]:newtree_end_idxs[newt]];
    array[N_max_newyears] real spring_temps_tree = newspring_temps[newtree_start_idxs[newt]:newtree_end_idxs[newt]];
    
    array[N_max_newyears] vector[2] alpha;
    
    for (n in 1:N_max_newyears) {
      
      int global_idx = n + (N_max_newyears)*(newt-1);
      real tau_nm_m = inv_logit(logit(tau_nm_m0) + beta_nm_m * (prevsummer_temps_tree[n]-summertemp0));
      real tau_m_m = tau_m_m0;
      matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                           [1 - tau_m_m, tau_m_m] ];

      if(n == 1){
        alpha[n] = omega[,n] .* rho;
      }
      else{
        alpha[n] = omega[,n] .* (Gamma' * alpha[n - 1]);
      }
      alpha[n] /= max(alpha[n]);
      
      vector[2] lambda = alpha[n] / sum(alpha[n]);
      states_new[global_idx] = categorical_rng(lambda);
      
      states_new_plot[newtree_stand_idxs[newt], n] += states_new[global_idx]-1;
      
      if(states_new[global_idx] == 1) {
        if(bernoulli_rng(theta1)){
          seed_counts_new[global_idx] = 0;
        }else{
          seed_counts_new[global_idx] = neg_binomial_2_alt_log_safe_rng(log(lambda1), psi1);  //neg_binomial_alt_rng(lambda1, psi1)
        }
      }
      else if(states_new[global_idx] == 2) {
        real loglambda2 = log(lambda20) 
        + beta_lambda2_frost * (gdd_lastfrost_tree[n]-frostgdd0) 
        + beta_lambda2_spring * (spring_temps_tree[n]-springtemp0);
      
        seed_counts_new[global_idx] = neg_binomial_2_alt_log_safe_rng(loglambda2, psi2); // neg_binomial_alt_rng(lambda2, psi2) 
      }
      
      seed_counts_new_plot[newtree_stand_idxs[newt], n] += seed_counts_new[global_idx];
    }
  }
}
