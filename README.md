# BinaryReplicates

<!-- badges: start -->
[![R-CMD-check](https://github.com/pierrepudlo/BinaryReplicates/workflows/R-CMD-check/badge.svg)](https://github.com/pierrepudlo/BinaryReplicates/actions)
<!-- badges: end -->

Statistical methods for analyzing **binary replicates**, i.e., noisy binary measurements of latent binary states. This package implements the methods described in:

> Royer-Carenzi, M., Lorenzo, H., & Pudlo, P. (in press). Reconciling Binary Replicates: Beyond the Average. *Statistics in Medicine*.

## Overview

The package provides scoring functions to estimate the probability that an individual is in the positive state, given noisy replicated measurements:

| Method | Function | Requirements |
|--------|----------|--------------|
| Average-based | `average_scoring()` | None |
| Median-based | `median_scoring()` | None |
| MAP (EM algorithm) | `MAP_scoring()` | Fitted EM model |
| Likelihood-based | `likelihood_scoring()` | Known parameters |
| Bayesian | `bayesian_scoring()` | Fitted Bayesian model |

Additional features:
- **Classification** with inconclusive decisions (`classify_with_scores()`)
- **Prevalence estimation** from scores or Bayesian posterior
- **Credible intervals** for model parameters
- **Cross-validation** for model assessment (`cvEM()`)

## Statistical Model

For each individual $i$, we observe $n_i$ binary replicates. These are noisy measurements of a true latent state $T_i \in \{0, 1\}$:

$$
T_i \mid \theta \sim \text{Bernoulli}(\theta)
$$

$$
S_i \mid T_i, p, q \sim \text{Binomial}\big(n_i,\; T_i(1-q) + (1-T_i)p\big)
$$

where:
- $\theta \in (0, 1)$ is the **prevalence** (probability that $T_i = 1$)
- $p \in (0, 1/2)$ is the **false positive rate** (probability of observing 1 when $T_i = 0$)
- $q \in (0, 1/2)$ is the **false negative rate** (probability of observing 0 when $T_i = 1$)

The goal is to estimate the probability $\mathbb{P}(T_i = 1 \mid S_i = s_i)$ for each individual, which is given by:

$$
\mathbb{P}(T_i = 1 \mid S_i = s_i) = \frac{\theta \cdot (1-q)^{s_i} q^{n_i - s_i}}{\theta \cdot (1-q)^{s_i} q^{n_i - s_i} + (1-\theta) \cdot p^{s_i} (1-p)^{n_i - s_i}}
$$

## Installation

### Prerequisites

The package depends on `rstan` for Bayesian inference. Install it first by following the guide at:
https://mc-stan.org/install/

### Install from GitHub

```r
# install.packages("devtools")
devtools::install_github("pierrepudlo/BinaryReplicates")
```

## Quick Start

```r
library(BinaryReplicates)

# Load example data
data("periodontal")
ni <- periodontal$ni
si <- periodontal$si

# --- Fast approach: EM algorithm ---
fit_em <- EMFit(ni, si)
fit_em$parameters_hat
scores_MAP <- MAP_scoring(ni, si, fit_em)

# --- Full Bayesian approach ---
fit_bayes <- BayesianFit(ni, si, chains = 4, iter = 5000)
scores_Bayes <- bayesian_scoring(ni, si, fit_bayes)

# Classify individuals (0.5 = inconclusive)
classes_MAP <- classify_with_scores(scores_MAP, vL = 0.4, vU = 0.6)
classes_Bayes <- classify_with_scores(scores_Bayes, vL = 0.4, vU = 0.6)

# Compare classifications
table(MAP = classes_MAP, Bayesian = classes_Bayes)
```

## Usage Examples

### Generate synthetic data

```r
theta <- 0.4
p <- q <- 0.22
n <- 50
ni <- sample(2:6, n, replace = TRUE)
ti <- rbinom(n, 1, theta)
si <- rbinom(n, ni, ti * (1 - q) + (1 - ti) * p)
```

### Simple scoring methods

These methods require no model fitting:

```r
# Average-based scores
Y_A <- average_scoring(ni, si)
theta_A <- prevalence_estimate(Y_A)

# Median-based scores
Y_M <- median_scoring(ni, si)
theta_M <- prevalence_estimate(Y_M)
```

### MAP estimation with the EM algorithm

The EM algorithm estimates model parameters without full Bayesian inference:

```r
fit_em <- EMFit(ni, si)
fit_em$parameters_hat

# MAP scores use the estimated parameters
Y_MAP <- MAP_scoring(ni, si, fit_em)
theta_MAP <- prevalence_estimate(Y_MAP)

# Classification with thresholds
T_MAP <- classify_with_scores(Y_MAP, vL = 0.4, vU = 0.6)
```

### Full Bayesian inference

For full posterior inference, use Stan via `BayesianFit()`:
```r
# Fit the Bayesian model (uses Stan MCMC)
fit <- BayesianFit(ni, si, chains = 4, iter = 5000)
print(fit, pars = c("theta", "p", "q"))

# Credible intervals
credint(fit, level = 0.90)

# Bayesian scores and prevalence
Y_B <- bayesian_scoring(ni, si, fit)
theta_B <- bayesian_prevalence_estimate(fit)

# Classification
T_B <- classify_with_scores(Y_B, vL = 0.4, vU = 0.6)
```

To use multiple CPU cores for faster sampling:
```r
options(mc.cores = parallel::detectCores())
```

### Likelihood-based scores (known parameters)

When the true parameters are known (e.g., in simulations):

```r
Y_L <- likelihood_scoring(ni, si, list(theta = theta, p = p, q = q))
T_L <- classify_with_scores(Y_L, vL = 0.4, vU = 0.6)
```

### Prediction on new data

Compute predictive Bayesian scores for new observations:

```r
new_ni <- rep(10, 11)
new_si <- 0:10
new_scores <- predict_scores(new_ni, new_si, fit)
```

### Cross-validation

Assess model performance with cross-validation:

```r
cv_result <- cvEM(ni, si, N_cv = 10)
cv_classified <- classify_with_scores_cvEM(cv_result, ti = ti, vL = 0.4)
cv_classified$risk  # Empirical risk
```

## Documentation

For more details, see the package vignette:

```r
vignette("introduction", package = "BinaryReplicates")
```

## Citation

If you use this package, please cite:

```
Royer-Carenzi, M., Lorenzo, H., & Pudlo, P. (in press).
Reconciling Binary Replicates: Beyond the Average.
Statistics in Medicine.
```

## License
GPL (>= 3)