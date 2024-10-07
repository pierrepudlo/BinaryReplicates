#' Compute predictive Bayesian scores
#'
#' @param newdata_ni Numeric vector of the total numbers of replicates per individuals
#' @param newdata_si Numeric vector of the numbers of positive replicates per individuals
#' @param fit The `stanfit` object return by [BayesianFit]
#'
#' @importFrom rstan extract
#' @importFrom dplyr summarise group_by ungroup mutate
#' @importFrom magrittr %>%
#' @return
#' The `predict_scores` function returns the predictive Bayesian scores in a numeric vector.
#' The predictive Bayesian scores are the posterior probabilities that the true
#' latent \eqn{T_i}'s are equal to 1 on new data, averaged over the posterior distribution.
#' @export
#' @examples
#' data("periodontal")
#' theta <- mean(periodontal$ti)
#' fit <- BayesianFit(periodontal$ni, periodontal$si, chains = 2, iter = 5000)
#' credint(fit)
#' newdata <- data.frame(ni = rep(4, 5), si = 0:4)
#' newdata$Y_BP <- predict_scores(newdata$ni, newdata$si, fit)
#' newdata
#'
#' ## EM comparaison
#' fitEM <- EMFit(periodontal$si,periodontal$ni)
#' newdata$Y_EM <- likelihood_scoring(
#'   newdata$ni, newdata$si,
#'   fitEM$parameters_hat$theta, fitEM$parameters_hat$p, fitEM$parameters_hat$q)
predict_scores <- function(newdata_ni, newdata_si, fit) {
  if (length(newdata_ni) != length(newdata_si)) {
    stop("newdata_ni and newdata_ni must have the same length")
  }
  n <- length(newdata_ni)
  # Extract the latent variable
  post <- rstan::extract(fit, pars = c("theta", "p", "q"))
  n_post <- length(post$theta)
  # Compute the posterior probability
  newdata <- data.frame(
    id = rep(1:n, each = n_post),
    ni = rep(newdata_ni, each = n_post),
    si = rep(newdata_si, each = n_post),
    theta = rep(post$theta, times = n),
    p = rep(post$p, times = n),
    q = rep(post$q, times = n)
  ) %>%
    mutate(Y_B = likelihood_scoring(ni, si, list(theta=theta, p=p, q=q)))
  out <- newdata %>%
    group_by(id) %>%
    summarise(Y_B = mean(Y_B), .groups="keep")
  return(out$Y_B)
}
