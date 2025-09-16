// NOT useable code

data {
int<lower=1> C; // number of colonies 
int<lower=1> K; // number of traps
int<lower=1> T; // number of timepoints
int<lower=1> L; // number of landscapes
matrix[K, 2] trap_pos; // trap coordinates 
int y[C, K]; // counts of bees in traps
real lowerbound; // lower spatial bound on colony location 
real upper_y; // upper spatial bound on colony location (y)
real upper_x; // upper spatial bound on colony location (x)
matrix[K,T] floral; // floral quality at traps
}

parameters {
real<lower=0> rho; // foraging length scale
real theta; // floral resource attractiveness
array[C] real<lower=lowerbound, upper=upper_x> delta_x; // colony x coords
array[C] real<lower=lowerbound, upper=upper_y> delta_y; // colony y coords
vector[K] eps; // trap specific intercept
vector[C] zeta; // colony specific intercept
real<lower=0> sigma; // trap specific intercept variance
real<lower=0> tau; // colony specific intercept variance
}


transformed parameters {
  // for non-centered parametrization
  real<lower=0> tau_sqrt = sqrt(tau);
  real<lower=0> sigma_sqrt = sqrt(sigma); 
  vector[C] zeta_scale = zeta * tau_sqrt;
  vector[K] eps_scale = eps*sigma_sqrt;
}
  
  
model {
  
// temporary declarations
matrix[C,K] dis; // distance from colony C to trap K
matrix[C,K] lambda; // rate of captures for colony C at trap K
vector[K] multi_probs; // multinomial probability of capture from colony C at trap K

// priors
rho ~ lognormal(3.5, 0.5);
theta ~ normal(0, 1);
sigma ~ normal(0, 1);
tau ~ normal(0, 1); 
eps ~ normal(0, 1); 
zeta ~ normal(0, 1);

// calculate intensity
  for (l in 1:L) { // iterate over landscapes
    for (k in 1:K) { // iterate over traps
      for (i in 1:C) { // iterate over colonies
        dis[i, k] = sqrt(square(delta[i, 1] - trap[k,1]) + square(delta[i, 2] - trap[k,2]));
      
        for (t in 1:T){ // iterate over timepoints
          lambda[i, k] = (-dis[i,k]/rho) + theta*floral[k] + zeta_scale[i] + eps_scale[k];
          fq[trapnum, t] = normal_rng(0,1);
          yobs[i, trapnum, t] = poisson_log_rng(alpha - 0.5*d / square(rho) + theta*fq[trapnum, t]);
          // yobs[i, trapnum, t] = poisson_log_rng(alpha - 0.5*d / (square(rho)*exp(theta*fq[trapnum, t])));
        }
      }
    }
  }
}

