#' Counts and ERCC matrix for sorted cell lines.
#'
#' The HOS, HCT116, and A375 cell lines were FACS sorted, either as
#' singlets or multiplets, and sequenced. Multiplets were sequenced into known
#' combinations in order to know which connections are contained in each
#' multiplet. Note that many of the ERCC read fractions did not look good for
#' many of these cells and the sort and sequencing was repeated. The results
#' from the repeated experiment are in the countsSorted2 dataset.
#'
#' @docType data
#' @name countsSorted1
#' @format Matrix containing counts with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @return matrix
#' @examples
#' data(countsSorted1)
#'
NULL

#' Counts and ERCC matrix for sorted cell lines.
#'
#' The HOS, HCT116, and A375 cell lines were FACS sorted, either as
#' singlets or multiplets, and sequenced. Multiplets were sequenced into known
#' combinations in order to know which connections are contained in each
#' multiplet. Note that many of the ERCC read fractions did not look good for
#' many of these cells and the sort and sequencing was repeated. The results
#' from the repeated experiment are in the countsSorted2 dataset.
#'
#' @docType data
#' @name countsSortedERCC1
#' @format Matrix containing ERCC counts, with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @return matrix
#' @examples
#' data(countsSortedERCC1)
#'
NULL

