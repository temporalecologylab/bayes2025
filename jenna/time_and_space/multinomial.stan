//Modified from Pope & Jha, 2017, Cons. Genetics

data {
int<lower=1> C; // number of colonies 
int<lower=1> K; // number of traps
matrix[K, 2] trap_pos; // trap coordinates 
int y[C, K]; // counts of bees in traps
real lowerbound; // lower spatial bound on colony location 
real upper_y; // upper spatial bound on colony location (y)
real upper_x; // upper spatial bound on colony location (x)
vector[K] floral; // floral quality at traps
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
    for(i in 1:C){
      for(k in 1:K){
        dis[i, k] = sqrt(square(delta_x[i] - trap_pos[k,1]) + square(delta_y[i] - trap_pos[k,2]));
        lambda[i, k] = (-dis[i,k]/rho) + theta*floral[k] + zeta_scale[i] + eps_scale[k];
      }
      
      // compute  multinomial probabilities for colony i
      multi_probs = softmax(lambda[i,]');
    
      // add to target likelihood (trap counts for colony i)
      y[i,] ~ multinomial(multi_probs);
    }
}

