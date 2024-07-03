# BinaryReplicates

This is a package that implements methods described in the paper XXX.

As described in the paper, the package provides functions to compute:

- the average-based scorings and classifications
- the median-based scorings and classifications
- the likelihood-based scorings and classifications
- the bayesian-based scorings and classifications


## Dependencies

The package depends on `rstan`that can be installed following the guide written [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).


## Installation

You can install the package using the following command, after installing the `devtools` package:

```{r}
devtools::install_github("pierrepudlo/BinaryReplicates")
```

## Usage

The function to fit the Bayesian model is `BayesianFit`.
Here is an example on how to use it, based on data generated from the model
as follows.

```{r}
theta <- .4
p <- q <- .22
n <- 20
ni <- sample(2:6, n, replace = TRUE)
ti <- rbinom(n, 1, theta)
si <- rbinom(n, ni, ti*(1-q) + (1-ti)*p)
synth_data <- data.frame(ni = ni, si = si, ti=ti)
```

### Average- and median-based computations

The average-based scores, classifications and prevalence estimate can be computed with

```{r}
Y_A <- average_scoring(ni, si) # scoring
T_A <- classify_with_scores(Y_A, vL = .4, vU = .6) # classify
theta_A <- prevalence_estimate(Y_A) # prevalence estimate
```

Likewise for the median-based statistics with

```{r}
Y_M <- median_scoring(ni, si) # scoring
T_M <- classify_with_scores(Y_M, vL = .4, vU = .6) # classify
theta_M <- prevalence_estimate(Y_M) # prevalence estimate
```

### The likelihood-based scores, classification

The likelihood-based computations need to know values of the fixed parameters $\theta$, $p$ and $q$. It is useless to compute a likelihood-based prevalence estimate, as $\theta$ is supposed to be known. The likelihood-based scorings and classifications can be computed with

```{r}
Y_L <- likelihood_scoring(ni, si, theta, p, q) # scorings
T_L <- classify_with_scores(Y_L, vL = .4, vU = .6) # classifications
```

### Bayesian scorings, classifications and prevalence estimate

First, we need to fit the Bayesian model to the synthetic data with

```{r}
fit <- BayesianFit(ni, si, chains = 4, iter = 5000)
print(fit)
```

Note that if you have enough RAM and a multicore CPU on your machine, you can configure `rstan` to take advantage of them. Here is an example to use 4 cores :

```{r}
options(mc.cores = 4)
```


If you have installed the package `shinystan`, you can explore the posterior distribution with

```{r}
shinystan::launch_shinystan(fit)
```

Credible intervals of probability 80% can be obtained with

```{r}
credint(fit, level = .8)
```

The Bayesian scores, classifications and prevalence estimate can be obtained with

```{r}
Y_B <- bayesian_scoring(ni, si, fit) # scorings
T_B <- classify_with_scores(Y_B, vL = .4, vU = .6) # classifications
theta_B <- bayesian_prevalence_estimate(fit) # prevalence estimate
```

### Summary of the classifications

The following code provides a summary of the classifications obtained with the different methods. We count how many times each method classifies the data as $0$, $1/2$ or $1$, depending on the true value of $T$.

```{r}
library(tidyverse)
confusion <- synth_data %>%
  mutate(
    Status = ifelse(ti==1, "T=1", "T=0"),
    Averge = T_A, Median = T_M, Likelihood = T_L, Bayesian = T_B) %>%
  pivot_longer(cols = c(Averge, Median, Likelihood, Bayesian), 
               names_to = "Method", values_to = "Decision") %>%
  group_by(Status, Method, Decision) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(count = as.integer(count)) %>%
  pivot_wider(names_from = Decision, values_from = count, values_fill = 0) 
confusion
```
