## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BinaryReplicates)
library(tidyverse)

## -----------------------------------------------------------------------------
theta <- .4
p <- q <- .22
n <- 50
ni <- sample(2:6, n, replace = TRUE)
ti <- rbinom(n, 1, theta)
si <- rbinom(n, ni, ti*(1-q) + (1-ti)*p)
synth_data <- data.frame(ni = ni, si = si, ti=ti)

## -----------------------------------------------------------------------------
fit <- BayesianFit(ni, si, chains = 4, iter = 5000, refresh = 0)
print(fit, pars = c("theta", "p", "q"))

## -----------------------------------------------------------------------------
mle <- EMFit(ni, si)

## -----------------------------------------------------------------------------
synth_data <- synth_data %>% 
  mutate(Y_A = average_scoring(ni, si),
         Y_M = median_scoring(ni, si),
         Y_L = likelihood_scoring(ni, si, theta, p, q),
         Y_E = likelihood_scoring(ni, si, mle['theta'], mle['p'], mle['q']),
         Y_B = bayesian_scoring(ni, si, fit))

## -----------------------------------------------------------------------------
v_L <- .35
v_U <- 1 - v_L
synth_data <- synth_data %>% 
  mutate(T_A = classify_with_scores(Y_A, v_L, v_U),
         T_M = classify_with_scores(Y_M, v_L, v_U),
         T_L = classify_with_scores(Y_L, v_L, v_U),
         T_E = classify_with_scores(Y_E, v_L, v_U),
         T_B = classify_with_scores(Y_B, v_L, v_U))

## -----------------------------------------------------------------------------
confusion <- synth_data %>%
  mutate(
    Status = ifelse(ti==1, "T=1", "T=0"),
    Averge = T_A, Median = T_M, Likelihood = T_L, EM = T_E, Bayesian = T_B) %>%
  pivot_longer(cols = c(Averge, Median, Likelihood, EM, Bayesian), 
               names_to = "Method", values_to = "Decision") %>%
  group_by(Status, Method, Decision) %>%
  summarise(count = n(), .groups = "keep") %>%
  ungroup() %>%
  mutate(count = as.integer(count)) %>%
  pivot_wider(names_from = Decision, values_from = count, values_fill = 0) 
confusion

## -----------------------------------------------------------------------------
new_data <- data.frame(ni = rep(8, 9), si = 0:8)

## -----------------------------------------------------------------------------
new_data <- new_data %>% 
  mutate(Y_B = predict_scores(ni, si, fit),
         T_B = classify_with_scores(Y_B, v_L, v_U))
round(new_data, digits = 4)

## -----------------------------------------------------------------------------
new_data <- new_data %>% 
  mutate(Y_E = likelihood_scoring(ni, si, mle['theta'], mle['p'], mle['q']),
         T_E = classify_with_scores(Y_E, v_L, v_U))
round(new_data, digits = 4)

