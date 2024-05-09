# Save this file as `R/BayesianFit.R`

#' Fit the Bayesian model for Binary Replicates
#'
#' @export
#' @param ni Numeric or integer vector of n_i values.
#' @param si Numeric or integer vector of S_i values.
#' @param prior A list of prior parameters for the model.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
BayesianFit <- function(ni, si,
                    prior = list(a_FP=2, b_FP=2,
                                 a_FN=2, b_FN=2,
                                 a_T=.5, b_T=.5),...) {
  # Check the data
  if (!is.numeric(ni) | !is.numeric(si)) {
    stop("ni and si must be numeric vectors")
  }
  if (length(ni) != length(si)) {
    stop("ni and si must have the same length")
  }
  if (!is.integer(ni)) ni <- as.integer(ni)
  if (!is.integer(si)) si <- as.integer(si)

  n <- length(ni)
  standata <- list(n = n, ni = ni, si = si,
                   a_FP = prior$a_FP, b_FP = prior$b_FP,
                   a_FN = prior$a_FN, b_FN = prior$b_FN,
                   a_T = prior$a_T, b_T = prior$b_T)
  out <- rstan::sampling(stanmodels$BayesianModel, data = standata, ...)
  return(out)
}
