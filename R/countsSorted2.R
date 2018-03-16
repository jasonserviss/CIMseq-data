#' Counts, ERCC, and metadata for sorted cell lines.
#'
#' The HOS, HCT116, and A375 cell lines were FACS sorted, either as
#' singlets or multiplets, and sequenced (171116). Multiplets were sequenced
#' into known combinations in order to know which connections are contained in
#' each multiplet. This dataset is used in conjunction with the
#' inst/sortedMultiplets20171116/*.Rmd analyses in the sp.scRNAseqTesting
#' package.
#'
#' @title Counts data from sorted multiplets; experiment 2.
#' @docType data
#' @name countsSorted2
#' @format Two matrices, containing counts and ERCC counts, with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @keywords datasets
#' @examples
#' data(countsSorted2)
"countsSorted2"

#' Counts, ERCC, and metadata for sorted cell lines.
#'
#' The HOS, HCT116, and A375 cell lines were FACS sorted, either as
#' singlets or multiplets, and sequenced (171116). Multiplets were sequenced
#' into known combinations in order to know which connections are contained in
#' each multiplet. This dataset is used in conjunction with the
#' inst/sortedMultiplets20171116/*.Rmd analyses in the sp.scRNAseqTesting
#' package.
#'
#' @title ERCC counts data from sorted multiplets; experiment 2.
#' @docType data
#' @name countsSortedERCC2
#' @format Matrix containing ERCC counts, with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @keywords datasets
#' @examples
#' data(countsSortedERCC2)
#'
"countsSortedERCC2"

#' Counts, ERCC, and metadata for sorted cell lines.
#'
#' The HOS, HCT116, and A375 cell lines were FACS sorted, either as
#' singlets or multiplets, and sequenced (171116). Multiplets were sequenced
#' into known combinations in order to know which connections are contained in
#' each multiplet. This dataset is used in conjunction with the
#' inst/sortedMultiplets20171116/*.Rmd analyses in the sp.scRNAseqTesting
#' package.
#'
#' @title Counts data from sorted multiplets; experiment 2.
#' @docType data
#' @name countsSortedMeta2
#' @format Two matrices, containing counts and ERCC counts, with:
#' \describe{
#' \item{sample}{sample ID}
#' \item{plate}{plate ID}
#' \item{row}{Row position in plate}
#' \item{column}{Column position in plate}
#' \item{cellNumber}{Indicates if sample is a single cell or multiple cells}
#' \item{cellTypes}{Indicates samples cell types}
#' }
#' @keywords datasets
#' @examples
#' data(countsSortedMeta2)
"countsSortedMeta2"
