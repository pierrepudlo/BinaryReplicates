# BinaryReplicates

This is a package that implements methods described in the paper XXX.

## Dependencies

The package depends on `rstan`that can be installed following the guide written [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).


## Installation

You can install the package using the following command, after installing the `devtools` package:

```r
devtools::install_github("pierrepudlo/BinaryReplicates")
```

## Usage

The function to fit the Bayesian model is `BayesianFit`.
Here is an example on how to use it. 

First, let's generate some data according to the model:
```r
theta <- .4
p <- q <- .08
n <- 10
ni <- c(2, 3, 4, 2, 3, 4, 2, 3, 4, 2)
ti <- rbinom(n, 1, theta)
si <- rbinom(n, ni, ti*(1-q) + (1-ti)*p)
```

Then we can fit the model with the default non informative prior:

```r
library(BinaryReplicates)
fit <- BayesianFit(ni, si, chains = 1, iter = 5000)
print(fit)
```

If you have installed the package `shinystan`, you can explore the posterior distribution with

```r
shinystan::launch_shinystan(fit)
```

Credible intervals of probability 90% can be obtained with

```r
credint(fit, alpha = .1)
```

And posterior probabilities that $T_i=1$ can be obtained with

```r
postprobT(fit)
```
