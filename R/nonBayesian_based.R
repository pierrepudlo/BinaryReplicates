#' Compute the average-, median- and likelihood-based scores
#'
#' @title Non-Bayesian scoring methods
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param param A list with 3 entries:
#'   \code{theta} The probability that \eqn{T=1}, i.e., the prevalence,
#'   \code{p} The false positive rate,
#'   \code{q} The false negative rate.
#' @param fit The object returned by [EMFit] containing the results of the EM algorithm
#' @return A numeric vector of the scores
#'
#' @note For likelihood-based scores, the values of \eqn{\theta}, \eqn{p} and
#' \eqn{q} are required. Consequently, likelihood scoring is not directly applicable
#' in practice without parameter estimates.
#'
#' @examples
#' data("periodontal")
#' Y_A <- average_scoring(periodontal$ni, periodontal$si)
#' Y_M <- median_scoring(periodontal$ni, periodontal$si)
#' # In order to compute the likelihood-based scores, we need to know theta,
#' # p and q which can be estimated in this example as follows:
#' theta_hat <- mean(periodontal$ti)
#' cat("The prevalence in the data is ", theta_hat, "\n")
#' p_hat <- with(periodontal, sum(si[ti == 0]) / sum(ni[ti == 0]))
#' q_hat <- with(periodontal, 1 - sum(si[ti == 1]) / sum(ni[ti == 1]))
#' Y_L <- likelihood_scoring(periodontal$ni, periodontal$si,
#'                           list(theta = theta_hat, p = p_hat, q = q_hat))
#'
#' @name non_bayesian_scoring
#'
#' @rdname non_bayesian_scoring
#' @export
average_scoring <- function(ni, si) {
  if (!is.numeric(ni) || !is.numeric(si)) {
    stop("ni and si must be numeric vectors")
  }
  if (length(ni) != length(si)) {
    stop("ni and si must have the same length")
  }
  si / ni
}

#' @rdname non_bayesian_scoring
#'
#' @examples
#' data("periodontal")
#' Y_M <- median_scoring(periodontal$ni, periodontal$si)
#'
#' @export
median_scoring <- function(ni, si) {
  if (!is.numeric(ni) || !is.numeric(si)) {
    stop("ni and si must be numeric vectors")
  }
  if (length(ni) != length(si)) {
    stop("ni and si must have the same length")
  }
  # Use 2*si == ni for exact integer comparison (avoids float issues)
  ifelse(2 * si == ni, 0.5, ifelse(2 * si > ni, 1, 0))
}


#' @rdname non_bayesian_scoring
#' @importFrom stats dbinom
#' @export
likelihood_scoring <- function(ni, si, param) {
  if (!is.numeric(ni) || !is.numeric(si)) {
    stop("ni and si must be numeric vectors")
  }
  if (length(ni) != length(si)) {
    stop("ni and si must have the same length")
  }
  if (!is.list(param) || !all(c("theta", "p", "q") %in% names(param))) {
    stop("param must be a list with entries 'theta', 'p', and 'q'")
  }
  theta <- param$theta
  p <- param$p
  q <- param$q
  # Use log-space computation to avoid division by zero and underflow
  log_numer <- log(theta) + dbinom(si, ni, 1 - q, log = TRUE)
  log_denom_term2 <- log(1 - theta) + dbinom(si, ni, p, log = TRUE)
  # log_sum_exp trick for numerical stability
  log_denom <- log_numer + log1p(exp(log_denom_term2 - log_numer))
  exp(log_numer - log_denom)
}

#' @rdname non_bayesian_scoring
#' @export
#'
#' @examples
#' data("periodontal")
#' fit <- EMFit(periodontal$ni, periodontal$si)
#' Y_MAP <- MAP_scoring(periodontal$ni, periodontal$si, fit)
#'
#' @seealso [EMFit]
MAP_scoring <- function(ni, si, fit) {
  if (!is.list(fit) || !"parameters_hat" %in% names(fit)) {
    stop("fit must be a list with a 'parameters_hat' component")
  }
  likelihood_scoring(ni, si, fit$parameters_hat)
}

#' Classification based on a thresholding of the scores
#'
#' @param scores Numeric vector of the scores, computed with
#'               [average_scoring], [median_scoring], [MAP_scoring] or [bayesian_scoring]
#' @param vL The lower threshold
#' @param vU The upper threshold
#' @return A numeric vector of the classification (where 0.5 = inconclusive)
#'
#' @details Each decision \eqn{\hat{t}_i} is taken according to the following rule:
#' \deqn{
#'         \hat{t}_i = \begin{cases}
#'          0 & \text{if } y_i < v_L,\\
#'          1/2 & \text{if } v_L \leq y_i \leq v_U,\\
#'          1 & \text{if } y_i > v_U,
#'          \end{cases}
#' }
#' where \eqn{y_i} is the score for individual \eqn{i}.
#'
#' @export
classify_with_scores <- function(scores, vL, vU) {
  if (!is.numeric(scores)) {
    stop("scores must be a numeric vector")
  }
  ifelse(scores < vL, 0, ifelse(scores > vU, 1, 0.5))
}

#' Compute the average-/median- or MAP-based prevalence estimates based on the scores
#'
#' @param scores Numeric vector of the scores, computed with
#'               [average_scoring], [median_scoring] or [MAP_scoring]
#' @return A numeric value of the prevalence estimate
#'
#' @note We have shown that the median-based prevalence estimator is better
#' than the average-based prevalence estimator in terms of bias, except when
#' the prevalence is in an interval \eqn{J}. The length of \eqn{J} is small when the
#' number of replicates is always large.
#'
#' @examples
#' data("periodontal")
#' theta <- mean(periodontal$ti)
#' Y_A <- average_scoring(periodontal$ni, periodontal$si)
#' Y_M <- median_scoring(periodontal$ni, periodontal$si)
#' fit <- EMFit(periodontal$ni, periodontal$si)
#' Y_MAP <- MAP_scoring(periodontal$ni, periodontal$si, fit)
#' hat_theta_A <- prevalence_estimate(Y_A)
#' hat_theta_M <- prevalence_estimate(Y_M)
#' hat_theta_MAP <- prevalence_estimate(Y_MAP)
#' cat("The average-based prevalence estimate is ", hat_theta_A, "\n")
#' cat("The median-based prevalence estimate is ", hat_theta_M, "\n")
#' cat("The MAP-based prevalence estimate is ", hat_theta_MAP, "\n")
#' cat("The prevalence in the dataset is ", theta, "\n")
#'
#' @export
prevalence_estimate <- function(scores) {
  if (!is.numeric(scores)) {
    stop("scores must be a numeric vector")
  }
  mean(scores)
}