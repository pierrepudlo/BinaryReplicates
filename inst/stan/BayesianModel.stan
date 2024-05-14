functions{
  // The predictive probability that T_i = 1, given the data and the parameters.
  real pred(int s, int n, real p, real q, real theta){
      real numerateur;
      real denominateur;
      numerateur = theta * exp(binomial_lpmf(s | n, 1-q));
      denominateur = numerateur +
                     (1-theta) * exp(binomial_lpmf(s | n, p));
    return numerateur / denominateur;
  }
}


// The data accepted by the model.
data {
  int n;
  int si[n];
  int ni[n];
  real a_FP;
  real b_FP;
  real a_FN;
  real b_FN;
  real a_T;
  real b_T;
}

// The parameters accepted by the model.
parameters {
  real<lower=0, upper=.5> p;
  real<lower=0, upper=.5> q;
  real<lower=0, upper=1> theta;
}

// The model to be estimated.
//
model {
  p ~ beta(a_FP, b_FP);
  q ~ beta(a_FN, b_FN);
  theta ~ beta(a_T, b_T);
  //target += 2*beta_lpdf(2*p | a_FP, b_FP);
  //target += 2*beta_lpdf(2*q | a_FN, b_FN);
  //target += beta_lpdf(theta | a_T, b_T);
  for(i in 1:n){
    target += log(theta * exp(binomial_lpmf(si[i] | ni[i], 1-q)) +
              (1-theta) * exp(binomial_lpmf(si[i] | ni[i], p)));
  }
}

// Output
generated quantities {
  int Ti[n];
  for(i in 1:n){
    Ti[i] = bernoulli_rng(pred(si[i], ni[i], p, q, theta));
  }
}

