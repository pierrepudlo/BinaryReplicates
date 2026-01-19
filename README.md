# BinaryReplicates

This package implements the methods described in the paper XXX.

It provides functions to compute:

- average-based scores and classifications
- median-based scores and classifications
- likelihood-based scores and classifications
- Bayesian scores, classifications, and prevalence estimates


## Dependencies

The package depends on `rstan`, which can be installed by following the guide at:
https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started


## Installation

After installing the `devtools` package you can install this package with:

```{r}
devtools::install_github("pierrepudlo/BinaryReplicates")
```

## Usage

The function to fit the Bayesian model is `BayesianFit`.
The following example uses synthetic data generated from the model. 
First, we generate the synthetic data:

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

Compute average-based scores, classifications, and prevalence estimates:

```{r}
Y_A <- average_scoring(ni, si) # scoring
T_A <- classify_with_scores(Y_A, vL = .4, vU = .6) # classify
theta_A <- prevalence_estimate(Y_A) # prevalence estimate
```

Compute median-based statistics similarly:

```{r}
Y_M <- median_scoring(ni, si) # scoring
T_M <- classify_with_scores(Y_M, vL = .4, vU = .6) # classify
theta_M <- prevalence_estimate(Y_M) # prevalence estimate
```

### Likelihood-based scores and classifications

Likelihood-based computations require known values for the fixed parameters $\theta$, $p$ and $q$. Computing a likelihood-based prevalence estimate is unnecessary because $\theta$ is assumed known. Compute likelihood-based scores and classifications with:

```{r}
Y_L <- likelihood_scoring(ni, si, theta, p, q) # scorings
T_L <- classify_with_scores(Y_L, vL = .4, vU = .6) # classifications
```

### Bayesian scores, classifications, and prevalence estimate

Fit the Bayesian model to the synthetic data:

```{r}
fit <- BayesianFit(ni, si, chains = 4, iter = 5000)
print(fit)
```

If you have sufficient RAM and multiple CPU cores, configure rstan to use them:

```{r}
options(mc.cores = 4)
```

If shinystan is installed, explore the posterior with:

```{r}
shinystan::launch_shinystan(fit)
```

Obtain 80% credible intervals with:

```{r}
credint(fit, level = .8)
```

Compute Bayesian scores, classifications and the prevalence estimate:

```{r}
Y_B <- bayesian_scoring(ni, si, fit) # scorings
T_B <- classify_with_scores(Y_B, vL = .4, vU = .6) # classifications
theta_B <- bayesian_prevalence_estimate(fit) # prevalence estimate
```

### Summary of classifications

The following summarizes how often each method classifies observations as $0$, $1/2$ or $1$, stratified by the true value of $T$.

```{r}
library(tidyverse)
confusion <- synth_data %>%
  mutate(
    Status = ifelse(ti==1, "T=1", "T=0"),
    Average = T_A, Median = T_M, Likelihood = T_L, Bayesian = T_B) %>%
  pivot_longer(cols = c(Average, Median, Likelihood, Bayesian), 
               names_to = "Method", values_to = "Decision") %>%
  group_by(Status, Method, Decision) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(count = as.integer(count)) %>%
  pivot_wider(names_from = Decision, values_from = count, values_fill = 0) 
confusion
```
