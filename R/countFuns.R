#' moveGenesToRownames
#'
#' Moves the first column of the counts data.frame to rownames and removes
#' the old column.
#'
#' @name moveGenesToRownames
#' @rdname moveGenesToRownames
#' @param counts data.frame; A data frame with counts data.
#' @return The counts data.frame is returned with the first column moved to
#' the rownames of the data.frame.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(LETTERS, a = runif(26, 1, 100))
#' moveGenesToRownames(counts)
#'
NULL
#' @export

moveGenesToRownames <- function(counts) {
  if(!length(dim(counts))) {
    stop("dim(counts) must have a positive length.")
  }
  if(class(counts) != "data.frame") {
    stop("Counts is not a data.frame")
  }
  if((class(counts[[1]]) != "character") & (class(counts[[1]]) != "factor")) {
    stop("The first column of counts is not class character or class factor.")
  }
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  return(counts)
}

#' convertCountsToMatrix
#'
#' Coerces the counts data.frame into a matrix.
#'
#' @name convertCountsToMatrix
#' @rdname convertCountsToMatrix
#' @param counts data.frame; A data frame with counts data.
#' @return The counts data.frame coerced into a matrix.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a = runif(26, 1, 100), b = runif(26, 1, 100))
#' convertCountsToMatrix(counts)
#'
NULL
#' @export

convertCountsToMatrix <- function(counts) {
  if(class(counts) != "data.frame") {
    stop("Counts is not a data.frame")
  }
  if(!all(sapply(counts, class) %in% c("numeric", "integer"))) {
    stop("Non-numeric columns detected. Should gene names be moved to rownames?")
  }
  as.matrix(counts)
}

#' labelSingletsAndMultiplets
#'
#' Adds the prefix "s." to colnames of columns that contain singlets and "m." to
#' columns that contain multiplets.
#'
#' @name labelSingletsAndMultiplets
#' @rdname labelSingletsAndMultiplets
#' @param counts data.frame; A data frame with counts data.
#' @param ids character; a character vector of regex statments used to indicate
#' colnames of columns that contain SINGLETS.
#' @return The counts data.frame with the appropriate prefix attatched to the
#'  column names. Throws a warning if no singlets are detected.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a = runif(26), b = runif(26), c = runif(26))
#' labelSingletsAndMultiplets(counts, ids = c("a", "b"))
#'
NULL
#' @export

labelSingletsAndMultiplets <- function(counts, ids) {
  if(is.null(colnames(counts))) {
    stop("is.null(colnames(counts)) returned TRUE.")
  }
  if(!all(sapply(counts, class) %in% c("numeric", "integer"))) {
    stop("Non-numeric columns detected. Should gene names be moved to rownames?")
  }
  if(class(ids) != "character") {
    stop("The id argument must be a character vector.")
  }
  
  bool <- sapply(ids, function(x) grepl(x, colnames(counts)))
  
  if(sum(bool) == 0) {
    warning("No singlets were detected. Is your id argument valid?")
  }
  
  colnames(counts) <- ifelse(
    rowSums(bool) > 0,
    paste("s.", colnames(counts), sep = ""),
    paste("m.", colnames(counts), sep = "")
  )
  return(counts)
}

#' removeHTSEQsuffix
#'
#' HTSeq adds the suffix ".htseq" to column names when it reports counts. This
#' function removes that suffix from the column names of the supplied counts
#' data.frame.
#'
#' @name removeHTSEQsuffix
#' @rdname removeHTSEQsuffix
#' @param counts data.frame; A data frame with counts data.
#' @return The counts data.frame with the ".htseq" suffix removed from the
#'  column names.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a.htseq = runif(26), b.htseq = runif(26))
#' removeHTSEQsuffix(counts)
#'
NULL
#' @export

removeHTSEQsuffix <- function(counts) {
  if(is.null(colnames(counts))) {
    stop("is.null(colnames(counts)) returned TRUE.")
  }
  if(all(colnames(counts) == paste0("V", 1:ncol(counts)))) {
    warning("Your colnames are V1..Vi. Are these sample names?" )
  }
  if(!any(grepl("\\.htseq", colnames(counts)))) {
    warning("Could not find the .htseq suffix in colnames(counts)")
  }
  colnames(counts) <- gsub("(.*)\\.htseq$", "\\1", colnames(counts))
  return(counts)
}

#' detectERCCreads
#'
#' Detects the 92 ERCC reads in the rownames of the counts data.frame. ERCC
#' reads must be named with the standard naming convention that matches the
#' regex "^ERCC\\-[0-9]*$".
#'
#' @name detectERCCreads
#' @rdname detectERCCreads
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row contains an ERCC read. A warning is issued if all 92
#' ERCC reads are not detected.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(runif(100), row.names = c(1:8, paste0("ERCC-", 1:92)))
#' detectERCCreads(counts)
#'
NULL
#' @export

detectERCCreads <- function(counts) {
  if(class(rownames(counts)) != "character") {
    stop("rownames(counts) is not of class character.")
  }
  if(all(rownames(counts) == as.character(1:nrow(counts)))) {
    m <- "rownames(counts) = 1:nrow(counts). Are gene names in rownames counts?"
    stop(m)
  }
  ercc <- grepl("^ERCC\\-[0-9]*$", rownames(counts))
  
  if(sum(ercc) != 92) {
    warning("Couldn't detect all ERCC reads.")
  }
  
  return(ercc)
}

#' detectNonGenes
#'
#' HTSeq outputs several "non-genes" in the counts output. These include:
#' "__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned",
#' "__alignment_not_unique". Since these are measurements of the quantification
#' quality and are often undesirable in downstream analysis, this function
#' detects them in the submitted counts data.frame to allow easy removal.
#'
#' @name detectNonGenes
#' @rdname detectNonGenes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row contains a non-gene.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(runif(10), row.names = c(1:9, "__no_feature"))
#' detectNonGenes(counts)
#'
NULL
#' @export

detectNonGenes <- function(counts) {
  if(class(rownames(counts)) != "character") {
    stop("rownames(counts) is not of class character.")
  }
  if(all(rownames(counts) == as.character(1:nrow(counts)))) {
    m <- "rownames(counts) = 1:nrow(counts). Are gene names in rownames counts?"
    stop(m)
  }
  
  nonGenes <- c(
    "__no_feature", "__ambiguous", "__too_low_aQual",
    "__not_aligned", "__alignment_not_unique"
  )
  rownames(counts) %in% nonGenes
}

#' detectLowQualityGenes
#'
#' In gene expression counts data if can often be the case that some genes are
#' not detected. This can simply be due to the fact that the gene is not
#' expressed in the tissue or, in addition, that the sequencing depth was not
#' sufficient to detect the gene. In addition, some genes may be detected but in
#' so few samples or at such a low level that it makes the quantified value
#' highly unreliable. In these cases, it is desireable to remove the gene before
#' downstream analysis which is facilitated by this function.
#'
#' @name detectLowQualityGenes
#' @rdname detectLowQualityGenes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @param mincount numeric; A minimum rowSum for which rows with a higher rowSum
#' will be detected. Default = 0.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row meets both tested conditions.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(c(0, runif(10)), c(0, runif(10)), c(0, runif(10)))
#' detectLowQualityGenes(counts)
#'
NULL
#' @export

detectLowQualityGenes <- function(
  counts,
  mincount = 0
){
  #checks
  rowSums(counts) > mincount
}

#' detectLowQualityCells
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality, in many cases due to issues during the sample preperation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this using two methods.
#' First, the mincount argument detects samples whose sum across all genes is >
#' mincount. Second, we utilize a house keeping gene and assume its expression
#' to be normally distributed. We then detect samples where the probability of
#' the expression for the house keeping gene in that sample is greater than the
#' quantile.cut argument.
#'
#' @name detectLowQualityCells
#' @rdname detectLowQualityCells
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param mincount numeric; A minimum colSum for which columns with a higher
#' colSum will be detected. Default = 4e5.
#' @param geneName character; The gene name to use for the quantile cutoff. This
#' must be present in the rownames of the counts argument. Default is ACTB.
#' @param quantile.cut numeric; This indicates probability at which the quantile
#' cutoff will be calculated using the normal distribution.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with colSums > mincount.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(
#'  runif(2e4),
#'  runif(2e4, 1, 100),
#'  row.names = paste0(letters, 1:2e4)
#' )
#' detectLowQualityCells(counts, geneName = "a1")
#'
NULL
#' @export

detectLowQualityCells <- function(
  counts,
  mincount = 4e5,
  geneName = 'ACTB',
  quantileCut = 0.001
){
  #input checks
  
  #colsums check
  cs <- colSums(counts) > mincount
  
  #house keeping check
  counts.log <- .norm.log.counts(counts)
  cl.act <- counts.log[geneName,]
  cl.act.m <- median(cl.act)
  cl.act.sd <- sqrt(
    sum((cl.act[cl.act > cl.act.m] - cl.act.m) ^ 2) /
    (sum(cl.act > cl.act.m) - 1)
  )
  my.cut <- qnorm(p = quantileCut, mean = cl.act.m, sd = cl.act.sd)
  
  cs & counts.log[geneName, ] > my.cut
}

#calcualtes log cpm
.norm.log.counts <- function(counts) {
  norm.fact <- colSums(counts)
  counts.norm <- t(apply(counts, 1, .norm, n = norm.fact))
  counts.log <- log2(counts.norm)
}

#calculates cpm on one row
.norm <- function(x, n) {
  x / n * 1000000 + 1
}
