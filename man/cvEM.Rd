% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMFit.R
\name{cvEM}
\alias{cvEM}
\title{Function to perform cross-validation in the EM algorithm}
\usage{
cvEM(
  ni,
  si,
  ti = NULL,
  N_cv = NULL,
  N_init = 20,
  maxIter = 1000,
  errorMin = 1e-07,
  prior = list(a_FP = 1, b_FP = 1, a_FN = 1, b_FN = 1)
)
}
\arguments{
\item{ni}{Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual}

\item{si}{Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual}

\item{ti}{Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual
If \code{NULL}, the EM algorithm is used to estimate the parameters. Default to \code{NULL}. See details}

\item{N_cv}{The number of folds. Default to 20}

\item{N_init}{The number of initializations if \code{ti} is not provided. Default to 20 if  observations otherwise default to the number of observations, corresponding to leave-one-out cross-validation.}

\item{maxIter}{The maximum number of iterations if EM algorithm is used. Default to 1e3}

\item{errorMin}{The minimum error computed if EM algorithm is used. Default to 1e-7}

\item{prior}{A list of prior parameters for the model. The prior distribution is as follows:}
}
\value{
A list with the following components:
\describe{
\item{models}{ A list of the models for each fold}
\item{predictions}{ A list of the predictions for each fold}
}
}
\description{
Function to perform cross-validation in the EM algorithm
}
\details{
This function chooses its algorithm according to what is provided in the \code{ti} argument:
\describe{
\item{\code{ti} is fully provided}{The function computes the \emph{Maximum-A-Posteriori}
estimate, with an explicit formula.}
\item{\code{ti} is not provided}{The function uses the EM algorithm
to estimate the parameters.}
\item{\code{ti} is partially provided}{The function uses the EM algorithm to
estimate the parameters.}
}
}
\examples{
data("periodontal")
modelCV <- cvEM(periodontal$ni,periodontal$si)

}
\seealso{
\link{classify_with_scores}, \link{EMFit}
}
