#' A periodontal dataset
#'
#'
#' Data from enzymatic diagnostic tests to detect two organisms
#' Treponema denticola and Bacteroides gingivalis.
#'
#' @docType data
#'
#' @usage data(periodontal)
#'
#' @format A data frame with 50 rows and 3 variables:
#' \describe{
#'  \item{ni}{Total numbers of sites tested for each individual}
#'  \item{si}{Number of positive tests for each individual}
#'  \item{ti}{Status of the individual (1=infected, 0=non-infected)}
#' }
#'
#' @keywords datasets
#'
#' @references Hujoel PP, Moulton LH, Loesche WJ. Estimation of sensitivity and
#' specificity of site-specific diagnostic tests. J Periodontal Res.
#' 1990 Jul;25(4):193-6. doi: 10.1111/j.1600-0765.1990.tb00903.x.
#' Erratum in: J Periodontal Res 1990 Nov;25(6):377. PMID: 2197400.
#' Moore et al. (2013) Genetics 195:1077-1086
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/2197400}{PubMed})
#'
#'
#' @examples
#' data(periodontal)
#' hat_prevalence <- mean(periodontal$si/periodontal$ni)
#' hat_prevalence
#' # should be compared to:
#' mean(periodontal$ti)
"periodontal"
