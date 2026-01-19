functions {
  // The predictive probability that T_i = 1, given the data and the parameters.
  real pred(int s, int n, real p, real q, real theta) {
    real log_num = log(theta) + binomial_lpmf(s | n, 1 - q);
    real log_denom = log_sum_exp(
      log_num,
      log1m(theta) + binomial_lpmf(s | n, p)
    );
    return exp(log_num - log_denom);
  }
}

// The data accepted by the model.
data {
  int<lower=1> n;
  array[n] int<lower=0> si;
  array[n] int<lower=1> ni;
  real<lower=0> a_FP;
  real<lower=0> b_FP;
  real<lower=0> a_FN;
  real<lower=0> b_FN;
  real<lower=0> a_T;
  real<lower=0> b_T;
}

// The parameters accepted by the model.
parameters {
  real<lower=0, upper=0.5> p;
  real<lower=0, upper=0.5> q;
  real<lower=0, upper=1> theta;
}

// The model to be estimated.
model {
  p ~ beta(a_FP, b_FP) T[0, 0.5];
  q ~ beta(a_FN, b_FN) T[0, 0.5];
  theta ~ beta(a_T, b_T);

  for (i in 1:n) {
    target += log_sum_exp(
      log(theta) + binomial_lpmf(si[i] | ni[i], 1 - q),
      log1m(theta) + binomial_lpmf(si[i] | ni[i], p)
    );
  }
}

// Output
generated quantities {
  array[n] int Ti;
  array[n] real posterior_prob;

  for (i in 1:n) {
    posterior_prob[i] = pred(si[i], ni[i], p, q, theta);
    Ti[i] = bernoulli_rng(posterior_prob[i]);
  }
}