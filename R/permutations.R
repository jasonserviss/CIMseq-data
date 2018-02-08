#' Permutations of the countsSorted2 dataset.
#'
#' The countsSorted2 dataset was used to run the sp.scRNAseq pipeline. The cell
#' types, which are known, see \code{help(countsSorted2)}, were renamed in the
#' spUnsupervised object. The \code{permuteSwarm} was then run and both the
#' spSwarm object and the output from the \code{permuteSwarm} function were
#' saved.
#'
#' @title spSwarm object with "real" data from permutation analysis.
#' @docType data
#' @name sObjPermutations
#' @format An spSwarm object containing the "real" results.
#' @usage data(sObjPermutations)
#' @keywords datasets
#' @examples
#' data(sObjPermutations)
"sObjPermutations"

#' Permutations of the countsSorted2 dataset.
#'
#' The countsSorted2 dataset was used to run the sp.scRNAseq pipeline. The cell
#' types, which are known, see \code{help(countsSorted2)}, were renamed in the
#' spUnsupervised object. The \code{permuteSwarm} was then run and both the
#' spSwarm object and the output from the \code{permuteSwarm} function were
#' saved.
#'
#' @title Permutation results from sorted multiplets experiment 2.
#' @docType data
#' @name permutations
#' @format A list containing the permutation results. Can be processed into tidy
#' format with the \code{tidyPermutationData} function.
#' @usage data(permutations)
#' @keywords datasets
#' @examples
#' data(permutations)
"permutations"
