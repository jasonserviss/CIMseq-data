#' Permutations of the countsSorted2 dataset.
#'
#' The countsSorted2 dataset was used to run the sp.scRNAseq pipeline. The cell
#' types, which are known, see \code{help(countsSorted2)}, were renamed in the
#' spUnsupervised object. The \code{permuteSwarm} was then run and both the
#' spSwarm object and the output from the \code{permuteSwarm} function were
#' saved.
#'
#' @docType data
#' @name permutations
#' @format An spSwarm object containing the "real" results. A list containing
#' the permutation results. The later can be processed into tidy format with the
#' \code{tidyPermutationData} function.
#' @usage permutations or sObjPermutations
#' @return The spSwarm object containing the "real" data of a list containing
#' the permuted data.
#' @examples
#' #data(sObjPermutations)
#' #data(permutations)
NULL
