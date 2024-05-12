#' Compute the average-, median- and likelihood-based scores
#'
#' @title Non-Bayesian scoring methods
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param theta The probability of that \eqn{T=1}, i.e., the prevalence
#' @param p The false positivity rate
#' @param q The false negativity rate
#' @return A numeric vector of the scores
#'
#' @examples
#' data("periodontal")
#' Y_A <- average_scoring(periodontal$ni, periodontal$si)
#' Y_M <- median_scoring(periodontal$ni, periodontal$si)
#' theta <- mean(periodontal$ti)
#' cat("The prevalence in the data is ", theta, "\n")
#' p <- mean(Y_A[periodontal$ti == 0])
#' q <- mean(1-Y_A[periodontal$ti == 1])
#' Y_L <- likelihood_scoring(periodontal$ni, periodontal$si, theta, p, q)
#'
#' @name non_bayesian_scoring
#'
#' @rdname non_bayesian_scoring
#' @export
average_scoring <- function(ni, si) {
  si / ni
}

#' @rdname non_bayesian_scoring
#' @export
median_scoring <- function(ni, si) {
  ifelse(si == ni/2, .5, ifelse(si > ni/2, 1, 0))
}


#' @rdname non_bayesian_scoring
#' @export
likelihood_scoring <- function(ni, si, theta, p, q) {
  theta*dbinom(si, ni, 1-q) /
    (theta*dbinom(si, ni, 1-q) + (1-theta)*dbinom(si, ni, p))
}

#' Classification based on a thresholding of the scores
#'
#' @param scores Numeric vector of the scores, computed with
#'               [average_scoring], [median_scoring] or [bayesian_scoring]
#' @param vL The lower threshold
#' @param vU The upper threshold
#' @return A numeric vector of the classification (where .5 = inconclusive)
#' @export
classify_with_scores <- function(scores, vL, vU) {
  ifelse(scores < vL, 0, ifelse(scores > vU, 1, .5))
}

#' Compute the average- and median-based prevalence estimates based on the scores
#'
#' @param scores Numeric vector of the scores, computed with
#'               [average_scoring] or [median_scoring]
#' @return A numeric value of the prevalence estimate
#'
#' @examples
#' data("periodontal")
#' theta <- mean(periodontal$ti)
#' Y_A <- average_scoring(periodontal$ni, periodontal$si)
#' hat_theta_A <- prevalence_estimate(Y_A)
#' cat("The average-based prevalence estimate is ", hat_theta_A, "\n")
#' cat("The prevalence in the dataset is ", theta, "\n")
#'
#' @export
prevalence_estimate <- function(scores) {
  mean(scores)
}
