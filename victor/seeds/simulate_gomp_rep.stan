// code to simulate data, with the reparam. gompertz

functions {
  real gompertz(real t, real tau1, real tau50) {
    return exp(-log(2)*(log(2)/log(100))^((t-tau50)/(tau50-tau1)));
  }
}

data {
  int<lower=1> Ntimes; // number of discrete intervals 
  int<lower=0> d[Ntimes]; // day of obs., d at Ntimes = Tmax
  
  int<lower=0> Nseeds;

  real<lower=0, upper=1> pv; // probability of viability
  real<lower=0> tau1; 
  real<lower=tau1> tau50;
}

transformed parameters {
  vector<lower=0>[Ntimes+1] pg; //germination prob
  
  pg[2] = pv*(gompertz(d[1], tau1, tau50) - gompertz(0, tau1, tau50));
  for (t in 2:Ntimes) {
    pg[t+1] = pv*(gompertz(d[t], tau1, tau50) - gompertz(d[t-1], tau1, tau50));
  }

  pg[1] = (1-pv)+pv*(1-gompertz(d[Ntimes], tau1, tau50));
  
  pg = pg / sum(pg); //  for stability

}


generated quantities {
  int<lower=0> y_sim[Ntimes+1];

  y_sim = multinomial_rng(pg, Nseeds);
}


