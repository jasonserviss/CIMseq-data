#' Fetal pancreas counts and ERCC data.
#'
#' @title Fetal pancreas counts data.
#' @docType data
#' @name countsFp
#' @format matrix
#' \describe{
#'     \item{rownames}{Gene names. ERCC identified with "^ERCC-[0-9]*"}
#'     \item{colnames}{Samples. Singlets prefixed with "s" and multiplets "m".}
#' }
#' @keywords datasets
#' @examples
#' data(countsFp)
"countsFp"

#' Fetal pancreas counts and ERCC data.
#'
#' @title Fetal pancreas ERCC counts data.
#' @docType data
#' @name countsErccFp
#' @format matrix
#' \describe{
#'     \item{rownames}{Gene names. ERCC identified with "^ERCC-[0-9]*"}
#'     \item{colnames}{Samples. Singlets prefixed with "s" and multiplets "m".}
#' }
#' @keywords datasets
#' @examples
#' data(countsErccFp)
"countsErccFp"
