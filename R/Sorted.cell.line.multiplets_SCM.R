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
#' @name SCM.counts
#' @format Two matrices, containing counts and ERCC counts, with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @keywords datasets
#' @examples
#' data(SCM.counts)
"SCM.counts"

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
#' @name SCM.countsERCC
#' @format Matrix containing ERCC counts, with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @keywords datasets
#' @examples
#' data(SCM.countsERCC)
#'
"SCM.countsERCC"

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
#' @name SCM.Meta
#' @format Tibble
#' @keywords datasets
#' @examples
#' data(SCM.Meta)
"SCM.Meta"
