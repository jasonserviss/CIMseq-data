#' Fetal pancreas counts and ERCC data.
#'
#' @title Fetal pancreas counts and ERCC reads data.
#' @docType data
#' @name fetalPancreasCounts
#' @format matrix
#' \describe{
#'     \item{rownames}{Gene names. ERCC identified with "^ERCC-[0-9]*"}
#'     \item{colnames}{Samples. Singlets prefixed with "s" and multiplets "m".}
#' }
#' @usage countsFp or countsErccFp
#' @return Matrix of counts.
#' @examples
#' data(countsFp)
#' data(countsErccFp)
NULL
