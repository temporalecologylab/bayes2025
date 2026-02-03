functions {
  real gompertz(real t, real beta, real tau) {
    return exp(-log2() * exp( -(2 / log2()) 
    * beta * (t - tau)));
  }
}

data {
  int<lower=1> Ntimes; // number of discrete intervals 
  int<lower=0> d[Ntimes]; // day of obs., d at Ntimes = Tmax
  int<lower=0> y[Ntimes+1]; // first element is #seeds that never germinated
}

parameters {
  real<lower=0, upper=1> pv; // probability of viability
  real<lower=0> beta; // beta > 0? that's my expectation
  real<lower=0> tau;
}

transformed parameters {
  vector<lower=0>[Ntimes+1] pg; //germination prob
  
  pg[2] = pv*(gompertz(d[1], beta, tau) - gompertz(0, beta, tau));
  for (t in 2:Ntimes) {
    pg[t+1] = pv*(gompertz(d[t], beta, tau) - gompertz(d[t-1], beta, tau));
  }

  pg[1] = (1-pv)+pv*(1-gompertz(d[Ntimes], beta, tau));
  
  pg = pg / sum(pg); //  for stability

}

model {
  // priors, au pif
  pv ~ beta(4,2);
  tau ~ normal(10,3);
  beta ~ normal(0,0.5/2.57);

  // likelihood
  y ~ multinomial(pg);
}

generated quantities {
  int<lower=0> y_pred[Ntimes+1];

  y_pred = multinomial_rng(pg, sum(y));
}
