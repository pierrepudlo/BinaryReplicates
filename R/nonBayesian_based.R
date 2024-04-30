#' Compute the average-based scorings
#'
#' @export
#' @param ni A numeric vector of the total number of replicates for each individual
#' @param si A numeric vector of the number of replicates equal to 1 for each individual
#' @return A numeric vector of the average-based scorings
average_based_scorings <- function(ni, si) {
  si / ni
}

#' Compute the median-based scorings
#'
#' @export
#' @param ni A numeric vector of the total number of replicates for each individual
#' @param si A numeric vector of the number of replicates equal to 1 for each individual
#' @return A numeric vector of the median-based scorings
median_based_scorings <- function(ni, si) {
  ifelse(si == ni/2, .5, ifelse(si > ni/2, 1, 0))
}

#' Compute the most likely-based scorings
#' @export
#' @param ni A numeric vector of the total number of replicates for each individual
#' @param si A numeric vector of the number of replicates equal to 1 for each individual
#' @param theta The probability of that T=1
#' @param p The false positive rate
#' @param q The false negative rate
#' @return A numeric vector of the most likely-based scorings
most_likely_based_scorings <- function(ni, si, theta, p, q) {
  theta*dbinom(si, ni, 1-q) /
    (theta*dbinom(si, ni, 1-q) + (1-theta)*dbinom(si, ni, p))
}

#' Classification based on thresholding the scorings
#' @export
#' @param scorings A numeric vector of the scorings
#' @param vL The lower threshold
#' @param vU The upper threshold
#' @return A numeric vector of the classification (where .5 = inconclusive)
classification_from_scoring <- function(scorings, vL, vU) {
  ifelse(scorings < vL, 0, ifelse(scorings > vU, 1, .5))
}

#' Compute the average,median-based prevalence estimate based on the scorings
#' @export
#' @param scorings A numeric vector of the scorings
#' @return A numeric value of the prevalence estimator
prevalence_estimate <- function(scorings) {
  mean(scorings)
}
