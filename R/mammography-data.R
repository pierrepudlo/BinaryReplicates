#' A mammography dataset
#'
#' Data from a mammography screening program.
#' The dataset is not the original dataset, but an imputed dataset based on
#' the summary statistics available publicly.
#'
#' @docType data
#'
#' @keywords datasets
#' @source Beam C, Conant E, Sickles E. Association of volume and volume-independent factors with accuracy in screening mammogram interpretation. \emph{JNCI.} 2003;95:282-290.
#' @name mammography-datasets
NULL

#' @rdname mammography-datasets
#' @format The data frame \code{mammography} has 148 rows and 3 variables
#' \describe{
#'   \item{ti}{True latent state of the individual (1=positive, 0=negative)}
#'   \item{ni}{Number of mammography tests performed for each individual}
#'   \item{si}{Number of positive mammography tests for each individual}
#' }
#'
#' @usage data(mammography)
"mammography"

#' @rdname mammography-datasets
#' @format The matrix \code{observed} is a 110x148 array with the observed replicates, one column
#' for each individual. The rows correspond to the replicates, and the values are either 0 or 1.
#' @usage data(observed)
"observed"
