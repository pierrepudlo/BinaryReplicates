#' A mammography dataset
#'
#' Data from a mammography screening program.
#'
#' @docType data
#'
#' @usage data(mammography, observed)
#'
#' @format A data frame \emph{mammography} with 148 rows and 3 variables
#' \describe{
#'   \item{ti}{True latent state of the individual (1=positive, 0=negative)}
#'   \item{ni}{Number of mammography tests performed for each individual}
#'   \item{si}{Number of positive mammography tests for each individual}
#' }
#' and an 110x148 array \emph{observed} with the observed replicates, one column
#' for each individual.
#'
#' @keywords datasets
#' @name mammography-datasets
NULL

#' @rdname mammography-datasets
"mammography"

#' @rdname mammography-datasets
"observed"
