#' Compute the *Maximum-A-Posteriori* estimate with the EM algorithm
#'
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param ti Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual.
#'  If \code{NULL}, the EM algorithm is used to estimate the parameters. Defaults to \code{NULL}. See details.
#' @param N_init The number of initializations if \code{ti} is not provided. Defaults to 20.
#' @param maxIter The maximum number of iterations if the EM algorithm is used. Defaults to 1e3.
#' @param errorMin The minimum error for convergence if the EM algorithm is used. Defaults to 1e-7.
#' @param prior A list of prior parameters for the model. The prior distribution is as follows:
#'
#' \itemize{
#' \item The false positive rate: \eqn{p \sim \text{Beta}(a_{FP}, b_{FP})}
#' \item The false negative rate: \eqn{q \sim \text{Beta}(a_{FN}, b_{FN})}
#' }
#' @return A list with the following components:
#' \describe{
#' \item{score}{ The estimated values of the scores}
#' \item{parameters_hat}{ The estimated values of the parameters \eqn{\theta}, \eqn{p} and \eqn{q}}
#' }
#'
#' @details
#' This function chooses its algorithm according to what is provided in the \code{ti} argument:
#' \describe{
#'   \item{\code{ti} is fully provided}{The function computes the \emph{Maximum-A-Posteriori}
#'   estimate, with an explicit formula.}
#'   \item{\code{ti} is not provided}{The function uses the EM algorithm
#'   to estimate the parameters.}
#'   \item{\code{ti} is partially provided}{The function uses the EM algorithm to
#'   estimate the parameters.}
#' }
#'
#' @importFrom stats rbinom dbinom
#'
#' @examples
#' data("periodontal")
#' # Get ML estimate knowing the true values of the latent ti's
#' periodontal_ml <- EMFit(periodontal$ni, periodontal$si, periodontal$ti)
#' # Get MAP estimate without knowing the true values of the latent ti's
#' periodontal_EM <- EMFit(periodontal$ni, periodontal$si, ti = NULL)
#'
#' @seealso [classify_with_scores]
#' @importFrom stats rbeta runif rbinom sd
#' @export
EMFit <- function(ni, si, ti = NULL,
                  prior = list(a_FP = 2, b_FP = 2, a_FN = 2, b_FN = 2),
                  N_init = 20, maxIter = 1e3, errorMin = 1e-7) {

  # Input validation
  if (!is.numeric(ni) || !is.numeric(si)) {
    stop("ni and si must be numeric vectors")
  }
  if (length(ni) != length(si)) {
    stop("ni and si must have the same length")
  }

  # Helper function to compute MAP estimates of p and q
  update_pq_estimates <- function(ni, si, score, prior) {
    p_hat <- (sum(si * (1 - score)) + prior$a_FP - 1) /
      (sum(ni * (1 - score)) + prior$a_FP + prior$b_FP - 2)
    q_hat <- (sum((ni - si) * score) + prior$a_FN - 1) /
      (sum(ni * score) + prior$a_FN + prior$b_FN - 2)
    list(p = p_hat, q = q_hat)
  }

  # Compute log-likelihood
  get_log_likelihood <- function(ni, si, p, q, theta) {
    n <- length(ni)
    likelihoods <- theta * (1 - q)^si * q^(ni - si) +
      (1 - theta) * p^si * (1 - p)^(ni - si)
    sum(log(likelihoods))
  }

  # Core EM algorithm
  em_bin <- function(ni, si, ti = NULL, maxIter = 1e3, errorMin = 1e-7, prior) {
    n <- length(si)

    if (!is.null(ti) && !any(is.na(ti))) {
      # Closed-form solution when ti is fully known
      theta_hat <- mean(ti)
      if (sum(ti < 1/2) == 0) {
        p_hat <- 0
      } else {
        p_hat <- sum(si[ti < 1/2]) / sum(ni[ti < 1/2])
      }
      if (sum(ti >= 1/2) == 0) {
        q_hat <- 0
      } else {
        q_hat <- sum((ni - si)[ti >= 1/2]) / sum(ni[ti >= 1/2])
      }
      out <- list(
        score = ti,
        ti_clas = ti,
        parameters_hat = data.frame(theta = theta_hat, p = p_hat, q = q_hat)
      )
    } else {
      # EM algorithm
      if (is.null(ti)) {
        score <- rbeta(n, shape1 = si, shape2 = ni - si)
        id_na <- 1:n
      } else {
        id_na <- which(is.na(ti))
        score <- ti
        score[id_na] <- rbeta(length(id_na),
                               shape1 = si[id_na] + 1,
                               shape2 = ni[id_na] - si[id_na] + 1)
      }

      iter <- 1
      error <- 2 * errorMin
      theta_hat <- mean(score)

      # Random initialization for p and q (Beta(2,2) prior mean = 0.5)
      p_hat <- rbeta(1, 2, 2)
      q_hat <- rbeta(1, 2, 2)

      while (error > errorMin && iter < maxIter) {
        # E-step: update scores (log-space for numerical stability)
        score_0 <- score
        log_numer <- log(theta_hat) + (ni - si) * log(q_hat) + si * log(1 - q_hat)
        log_denom_r <- log(1 - theta_hat) + si * log(p_hat) +
                       (ni - si) * log(1 - p_hat)
        # log_sum_exp trick: log(a + b) = log_a + log(1 + exp(log_b - log_a))
        log_denom <- log_numer + log1p(exp(log_denom_r - log_numer))
        score[id_na] <- exp(log_numer[id_na] - log_denom[id_na])

        # M-step: update parameters
        theta_hat <- mean(score)
        pq <- update_pq_estimates(ni, si, score, prior)
        p_hat <- pq$p
        q_hat <- pq$q

        # Handle edge case
        if (is.na(q_hat)) q_hat <- runif(1, 0, 1/2)

        # Ensure identifiability (p, q < 0.5)
        if (p_hat > 1/2 && q_hat > 1/2) {
          p_hat <- 1 - p_hat
          q_hat <- 1 - q_hat
          theta_hat <- 1 - theta_hat
          score[id_na] <- 1 - score[id_na]
        }

        error <- mean(abs(score_0[id_na] - score[id_na]))
        iter <- iter + 1
      }

      out <- list(
        score = score,
        parameters_hat = data.frame(theta = theta_hat, p = p_hat, q = q_hat)
      )
    }
    out
  }

  # Main logic
  if (!is.null(ti) && !any(is.na(ti))) {
    out <- em_bin(ni, si, ti, prior = prior, errorMin = errorMin, maxIter = maxIter)
  } else {
    result <- matrix(NA, N_init, 4)
    colnames(result) <- c("theta", "p", "q", "LogPosterior")
    allRes <- list()

    for (i in 1:N_init) {
      allRes[[i]] <- em_bin(ni, si, prior = prior, errorMin = errorMin, maxIter = maxIter)
      parai <- allRes[[i]]$parameters_hat
      loglikeli <- get_log_likelihood(ni, si, p = parai$p, q = parai$q, theta = parai$theta)
      log_prior <- (prior$a_FP - 1) * log(parai$p) + (prior$b_FP - 1) * log(1 - parai$p) +
                   (prior$a_FN - 1) * log(parai$q) + (prior$b_FN - 1) * log(1 - parai$q)
      logPoster <- loglikeli + log_prior
      result[i, ] <- c(as.numeric(parai), logPoster)
    }

    i_opt <- which.max(result[, 4])
    out <- allRes[[i_opt]]
  }
  out
}

#' Cross-validation for the EM algorithm
#'
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param ti Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual.
#'  If \code{NULL}, the EM algorithm is used to estimate the parameters. Defaults to \code{NULL}. See details.
#' @param N_cv The number of folds. Defaults to 20.
#' @param N_init The number of initializations if \code{ti} is not provided. Defaults to 20.
#' @param maxIter The maximum number of iterations if the EM algorithm is used. Defaults to 1e3.
#' @param errorMin The minimum error for convergence if the EM algorithm is used. Defaults to 1e-7.
#' @param prior A list of prior parameters for the model.
#'
#' @details
#' This function chooses its algorithm according to what is provided in the \code{ti} argument:
#' \describe{
#'   \item{\code{ti} is fully provided}{The function computes the \emph{Maximum-A-Posteriori}
#'   estimate, with an explicit formula.}
#'   \item{\code{ti} is not provided}{The function uses the EM algorithm
#'   to estimate the parameters.}
#'   \item{\code{ti} is partially provided}{The function uses the EM algorithm to
#'   estimate the parameters.}
#' }
#'
#' @return A list with the following components:
#' \describe{
#' \item{models}{ A list of the models for each fold}
#' \item{predictions}{ A list of the predictions for each fold}
#' }
#'
#' @examples
#' data("periodontal")
#' modelCV <- cvEM(periodontal$ni, periodontal$si)
#'
#' @seealso [classify_with_scores], [EMFit]
#' @export
cvEM <- function(ni, si, ti = NULL, N_cv = NULL,
                 N_init = 20, maxIter = 1e3, errorMin = 1e-7,
                 prior = list(a_FP = 2, b_FP = 2, a_FN = 2, b_FN = 2)) {

  n <- length(ni)
  if (is.null(N_cv)) N_cv <- min(20, n)
  if (N_cv < 2) stop("Choose at least 2 folds")
  if (N_cv %% 1 != 0) stop("Select an integer as number of folds")
  if (N_cv > n) stop("Select at most n folds")

  id_cv <- unlist(lapply(1:N_cv, rep, floor(n / N_cv)))
  if (n %% N_cv > 0) {
    id_cv <- c(id_cv, 1:(n %% N_cv))
  }
  id_cv <- sample(id_cv, size = n, replace = FALSE)

  MODELS <- list()
  PREDICTIONS <- list()

  for (i in 1:N_cv) {
    id_test <- id_cv == i
    id_train <- !id_test
    ni_train <- ni[id_train]
    si_train <- si[id_train]
    ni_test <- ni[id_test]
    si_test <- si[id_test]

    if (!is.null(ti)) {
      ti_train <- ti[id_train]
      ti_test <- ti[id_test]
    } else {
      ti_train <- NULL
      ti_test <- NULL
    }

    MODELS[[i]] <- EMFit(
      ni = ni_train, si = si_train, ti = ti_train,
      N_init = N_init, maxIter = maxIter, errorMin = errorMin,
      prior = prior
    )

    PREDICTIONS[[i]] <- list()
    PREDICTIONS[[i]]$scores_predicted <- likelihood_scoring(
      ni_test, si_test, MODELS[[i]]$parameters_hat
    )
    names(PREDICTIONS[[i]]$scores_predicted) <- which(id_test)

    if (!is.null(ti)) {
      PREDICTIONS[[i]]$class_observed <- ti_test
      names(PREDICTIONS[[i]]$class_observed) <- which(id_test)
    } else {
      PREDICTIONS[[i]]$class_observed <- NULL
    }
  }

  out <- list(models = MODELS, predictions = PREDICTIONS)
  class(out) <- "cvEM"
  out
}

#' Perform classification on the scores for each fold of a cvEM object
#'
#' Note that if \eqn{t_i} is provided, the empirical risk is estimated with
#' \eqn{a=v_L}.
#'
#' @param object An object of class cvEM
#' @param ti Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual.
#'  If \code{NULL}, the risk is not computed. Defaults to \code{NULL}.
#' @param vL The lower threshold for classification. Defaults to 0.5.
#' @param vU The upper threshold for classification. Defaults to 0.5.
#'
#' @return A cvEM object with the following components:
#' \describe{
#' \item{predictions}{ A list of the predictions for each fold}
#' \item{risk}{ The empirical risk if \eqn{t_i} is provided}
#' }
#' @export
#'
#' @examples
#' data(periodontal)
#' modelCV <- cvEM(periodontal$ni, periodontal$si)
#' modelCV2 <- classify_with_scores_cvEM(modelCV, vL = 0.4)
classify_with_scores_cvEM <- function(object, ti = NULL, vL = 0.5, vU = 0.5) {

  compute_risk <- function(ti_hat, ti, a) {
    n <- length(ti)
    risk <- 0
    for (i in 1:n) {
      if (abs(ti[i] - ti_hat[i]) == 1) risk <- risk + 1
      if (abs(ti[i] - ti_hat[i]) == 1/2) risk <- risk + a
    }
    risk
  }

  if (!inherits(object, "cvEM")) {
    stop("The object must be of class cvEM")
  }

  N_cv <- length(object$predictions)
  if (!is.null(ti)) RISK <- 0

  for (i in 1:N_cv) {
    object$predictions[[i]]$classes_predicted <- classify_with_scores(
      object$predictions[[i]]$scores_predicted, vL = vL, vU = vU
    )
    id_i <- as.numeric(names(object$predictions[[i]]$scores_predicted))
    names(object$predictions[[i]]$classes_predicted) <- id_i

    if (!is.null(ti)) {
      RISK <- RISK + compute_risk(
        object$predictions[[i]]$classes_predicted, ti[id_i], a = vL
      )
    }
  }

  if (!is.null(ti)) object$risk <- RISK
  object
}