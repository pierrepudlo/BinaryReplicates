# BinaryReplicates

This is a package that implements methods described in the paper XXX.

As described in the paper, the package provides functions to compute:

- the average-based scorings and classifications
- the median-based scorings and classifications
- the most likely-based scorings and classifications
- the bayesian-based scorings and classifications


## Dependencies

The package depends on `rstan`that can be installed following the guide written [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).


## Installation

You can install the package using the following command, after installing the `devtools` package:

```r
devtools::install_github("pierrepudlo/BinaryReplicates")
```

## Usage

The function to fit the Bayesian model is `BayesianFit`.
Here is an example on how to use it, based on data generated from the model
as follows.

```r
theta <- .4
p <- q <- .08
n <- 10
ni <- c(2, 3, 4, 2, 3, 4, 2, 3, 4, 2)
ti <- rbinom(n, 1, theta)
si <- rbinom(n, ni, ti*(1-q) + (1-ti)*p)
synth_data <- data.frame(ni = ni, si = si, ti=ti)
```

## Average- and median-based computations

The average-based scorings, classifications and prevalence estimate can be computed with

```r
Y_A <- average_based_scorings(ni, si) # scorings
T_A <- classification_from_scoring(Y_A, vL = .4, vU = .6) # classifications
theta_A <- prevalence_estimate(Y_A) # prevalence estimate
```

Likewise for the median-based statistics with

```r
Y_M <- median_based_scorings(ni, si) # scorings
T_M <- classification_from_scoring(Y_M, vL = .4, vU = .6) # classifications
theta_M <- prevalence_estimate(Y_M) # prevalence estimate
```

## Bayesian scorings, classifications and prevalence estimate

First, we need to fit the Bayesian model to the synthetic data with

```r
fit <- BayesianFit(ni, si, chains = 4, iter = 5000)
print(fit)
```

Note that if you have enough RAM and a multicore CPU on your machine, you can configure `rstan` to take advantage of them. Here is an example to use 4 cores :

```r
options(mc.cores = 4)
```


If you have installed the package `shinystan`, you can explore the posterior distribution with

```r
shinystan::launch_shinystan(fit)
```

Credible intervals of probability 90% can be obtained with

```r
credint(fit, alpha = .1)
```

The Bayesian scorings, classifications and prevalence estimate can be obtained with

```r
Y_B <- bayesian_scorings(fit) # scorings
T_B <- classification_from_scoring(Y_B, vL = .4, vU = .6) # classifications
theta_B <- bayesian_prevalence_estimate(fit) # prevalence estimate
```
