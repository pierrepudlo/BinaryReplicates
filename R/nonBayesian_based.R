#' Compute the average, median and most likely-based scorings
#'
#' @title Non-Bayesian scorings
#' @param ni A numeric vector of the total number of replicates for each individual
#' @param si A numeric vector of the number of replicates equal to 1 for each individual
#' @param theta The probability of that T=1, i.e., the prevalence
#' @param p The false positive rate
#' @param q The false negative rate
#' @return A numeric vector of the scorings
#'
#' @examples
#' data("periodontal")
#' Y_A <- average_based_scorings(periodontal$ni, periodontal$si)
#' Y_M <- median_based_scorings(periodontal$ni, periodontal$si)
#' theta <- mean(periodontal$ti)
#' cat("The prevalence in the data is ", theta, "\n")
#' p <- mean(Y_A[periodontal$ti == 0])
#' q <- mean(1-Y_A[periodontal$ti == 1])
#' Y_L <- most_likely_based_scorings(periodontal$ni, periodontal$si, theta, p, q)
#'
#' @name non_bayesian_scorings
#'
#' @rdname non_bayesian_scorings
#' @export
average_based_scorings <- function(ni, si) {
  si / ni
}

#' @rdname non_bayesian_scorings
#' @export
median_based_scorings <- function(ni, si) {
  ifelse(si == ni/2, .5, ifelse(si > ni/2, 1, 0))
}


#' @rdname non_bayesian_scorings
#' @export
most_likely_based_scorings <- function(ni, si, theta, p, q) {
  theta*dbinom(si, ni, 1-q) /
    (theta*dbinom(si, ni, 1-q) + (1-theta)*dbinom(si, ni, p))
}

#' Classification based on thresholding the scorings
#'
#' @param scorings A numeric vector of the scorings
#' @param vL The lower threshold
#' @param vU The upper threshold
#' @return A numeric vector of the classification (where .5 = inconclusive)
#' @export
classification_from_scoring <- function(scorings, vL, vU) {
  ifelse(scorings < vL, 0, ifelse(scorings > vU, 1, .5))
}

#' Compute the average,median-based prevalence estimate based on the scorings
#'
#' @param scorings A numeric vector of the scorings
#' @return A numeric value of the prevalence estimator
#'
#' @examples
#' data("periodontal")
#' theta <- mean(periodontal$ti)
#' Y_A <- average_based_scorings(periodontal$ni, periodontal$si)
#' hat_theta_A <- prevalence_estimate(Y_A)
#' cat("The average-based prevalence estimate is ", hat_theta_A, "\n")
#' cat("The prevalence in the dataset is ", theta, "\n")
#'
#' @export
prevalence_estimate <- function(scorings) {
  mean(scorings)
}
