functions {
  vector uniform_point_rng(real lowerbound,
                           real upperbound_x,
                           real upperbound_y) {
    vector[2] x;
    
    x[1] = uniform_rng(lowerbound, upperbound_x);
    x[2] = uniform_rng(lowerbound, upperbound_y);
    
    return x;
  }
  
}

data {
  int<lower=0> C; // Number of colonies
  int<lower=0> K;      // Number of traps
  int<lower=1> L; // Number of landscapes
  int<lower=1> T; // Number of timepoints
  real alpha; // Intercept
  real theta; // Floral resource effect 
}

transformed data {
  int<lower=0> landscapesize = 1500; // Size of each landscape
  int<lower=0> trapgridsize = 300; // Size of each trapping grid
  real gridsize = sqrt(K); // Number of rows/columns in trap grid
  real stepsize = trapgridsize/(gridsize-1); // Distance between traps in grid
  
  real<lower=0> lowerbound = 0; // lower limit for colony locations
  real<lower=0> upperbound_x = 3*landscapesize;   // upper limit for colony locations (x direction)
  real<lower=0> upperbound_y = (L/3)*landscapesize;   // upper limit for colony locations (y direction)
  
  real<lower=0> rho = 100;
}

generated quantities {
  array[C] vector[2] delta;
  array[K*L] vector[2] trap_pos;
  array[K*L, T] real fq;
  
  array[C, K*L, T] int<lower=0> yobs;
  array[C, K*L] real<lower=0> sq_dist;
  
  for (i in 1:C)
    delta[i] = uniform_point_rng(lowerbound, upperbound_x, upperbound_y);
  
  for (l in 1:L) { // iterate over landscapes
    for (k in 1:K) { // iterate over traps
    
      //define landscape position
      int lancol = ((l-1) % 3);
      real lanrow = ceil(l/3.0) - 1;
      
      // define trap position
      int trapnum = k + (l-1)*K;
      int trapcol = (k-1) % to_int(gridsize);
      real traprow = ceil(k/gridsize) - 1;
      
      // get trap coordinates
      trap_pos[trapnum][1] = lancol*landscapesize + (landscapesize - trapgridsize)/2 + (trapcol)*stepsize;
      trap_pos[trapnum][2] = lanrow*landscapesize + (landscapesize - trapgridsize)/2 + (traprow)*stepsize;

      
      for (i in 1:C) { // iterate over colonies
        real d = squared_distance(trap_pos[trapnum], delta[i]);
        sq_dist[i, trapnum] = d;
        
        for (t in 1:T){ // iterate over timepoints
          fq[trapnum, t] = normal_rng(0,1);
          yobs[i, trapnum, t] = poisson_log_rng(alpha - 0.5*d / square(rho) + theta*fq[trapnum, t]);
          // yobs[i, trapnum, t] = poisson_log_rng(alpha - 0.5*d / (square(rho)*exp(theta*fq[trapnum, t])));
        }
      }
    }
    
  }
  
}
