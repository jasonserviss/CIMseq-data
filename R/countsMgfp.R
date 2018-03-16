#' Counts, ERCC, and metadata from the GFP-Lgr5 mouse small intestine and colon.
#'
#' @title Counts data from the GFP-Lgr5 mouse small intestine and colon.
#' @docType data
#' @name countsMgfp
#' @format Matrix containing counts with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @keywords datasets
#' @examples
#' data(countsMgfp)
#'
"countsMgfp"

#' Counts, ERCC, and metadata from the GFP-Lgr5 mouse small intestine and colon.
#'
#' @title ERCC counts data from the GFP-Lgr5 mouse small intestine and colon.
#' @docType data
#' @name countsMgfpERCC
#' @format Matrix containing ERCC counts with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @keywords datasets
#' @examples
#' data(countsMgfpERCC)
#'
"countsMgfpERCC"

#' Counts, ERCC, and metadata from the GFP-Lgr5 mouse small intestine and colon.
#'
#' @title Metadata from the GFP-Lgr5 mouse small intestine and colon.
#' @docType data
#' @name countsMgfpMeta
#' @format Tibble with:
#' \describe{
#' \item{sample}{Sample ID}
#' \item{plate}{Plate ID}
#' \item{row}{Row position in plate}
#' \item{column}{Column position in plate}
#' \item{GFP}{Samples GFP status}
#' \item{mouseID}{mouse ID}
#' \item{tissue}{tissue from which the sample was derived}
#' \item{cellNumber}{Indicates if sample is a single cell or multiple cells}
#' }
#' @keywords datasets
#' @examples
#' data(countsMgfpMeta)
#'
"countsMgfpMeta"
