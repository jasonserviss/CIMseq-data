#' getMetadata
#'
#' Utilizes the EngeMetadata package to download and format the metadata
#' associated with a specific project folder. Downloads all metadata/plate info
#' present in the annotation folder by default.
#'
#' @name getMetadata
#' @rdname getMetadata
#' @param projectName character; The name of the project. Must correspond with
#' the project folder name located at EngeLab/data.
#' @return The metadata tibble containing the information in the project
#' annotation folder.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom googledrive drive_ls
#' @importFrom dplyr "%>%" mutate if_else select everything
#' @importFrom purrr map_dfr
#' @importFrom EngeMetadata metadata

getMetadata <- function(projectName) {
  cellNumber <- prefix <- unique_key <- Well <- NULL
  annotationPath <- file.path('data', projectName, 'annotation')
  plates <- googledrive::drive_ls(annotationPath, trashed = FALSE)$name
  meta <- purrr::map_dfr(plates, function(p) {
    EngeMetadata::metadata(p, annotationPath)
  })

  if("cellNumber" %in% colnames(meta)) {
    meta <-  meta %>%
      dplyr::mutate(prefix = dplyr::if_else(cellNumber == "Singlet", "s.", "m.")) %>%
      dplyr::mutate(sample = paste0(prefix, unique_key, ".", Well)) %>%
      dplyr::select(-prefix) %>%
      dplyr::select(sample, dplyr::everything())
  }
  return(meta)
}

#' getCountsData
#'
#' Downloads the count data, located at EngeLab/project name/raw_data, for a
#' project and concatenates it into a data.frame. The column containing the gene
#' name, to be named HGN, is expected in all counts files.
#'
#' @name getCountsData
#' @rdname getCountsData
#' @param projectName character; The name of the project. Must correspond with
#' the project folder name located at EngeLab/data.
#' @return A data.frame with the unfiltered counts data.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom googledrive drive_ls drive_download
#' @importFrom purrr map2 map reduce
#' @importFrom dplyr full_join
#' @importFrom utils read.table

getCountsData <- function(projectName) {

  rawPath <- file.path('data', projectName, 'raw_counts')
  countsFiles <- googledrive::drive_ls(rawPath)$name
  fullPaths <- file.path(rawPath, countsFiles)
  localPaths <- file.path(tempdir(), countsFiles)

  trash <- purrr::map2(fullPaths, localPaths, function(fp, lp) {
    googledrive::drive_download(fp, lp, overwrite = TRUE)
  })
  rm(trash)

  #load counts
  loaded <- purrr::map(localPaths, utils::read.table, header = TRUE, sep = "\t")
  purrr::reduce(loaded, dplyr::full_join, by = "HGN")
}

#' filterCountsData
#'
#' This is a wrapper function for \code{\link{detectLowQualityGenes}},
#' \code{\link{detectLowQualityCells.totalCounts}},
#' \code{\link{detectLowQualityCells.ERCCfrac}},
#' \code{\link{detectLowQualityCells.housekeeping}},
#' \code{\link{convertCountsToMatrix}}. It uses the
#' "cellNumber" column in the metadata to annotate singlets and multiplets and,
#' therefore, the column must be present in the metadata. In addition, singlets
#' should be named "Singlet" in the metadata. The function also checks for the
#' presence of NAs in the counts data and returns an error if NA values are
#' detected.
#'
#' @name filterCountsData
#' @rdname filterCountsData
#' @param counts data.frame; The formated counts data.
#' @param countsERCC data.frame; The formated counts ERCC data.
#' @param filters Character; A vector indicating the types of filtering to be
#' performed. Options include: "genes", "totalCounts",
#' "ERCCfrac", and "housekeeping".
#' @param geneMinCount See \code{\link{detectLowQualityGenes}}
#' @param cellMinCount See \code{\link{detectLowQualityCells.totalCounts}}
#' @param geneName See \code{\link{detectLowQualityCells.housekeeping}}
#' @param quantileCut See \code{\link{detectLowQualityCells.housekeeping}}
#' @param percentile See \code{\link{detectLowQualityCells.ERCCfrac}}.
#' @return A list with the filtered counts as the first element and the filtered
#' ERCC genes as the second element.
#' @author Jason Serviss
#'
NULL
#' @export

filterCountsData <- function(
  counts, countsERCC,
  filters = c("genes", "totalCounts", "ERCCfrac", "housekeeping"),
  geneMinCount, cellMinCount, geneName, quantileCut.hk, quantileCut.ercc
){
  #remove low quality genes
  if("genes" %in% filters) {
    counts <- counts[detectLowQualityGenes(counts, mincount = geneMinCount), ]
  }

  #remove low quality cells
  lqc.totalCounts <- lqc.housekeeping <- lqc.ERCCfrac <- rep(TRUE, ncol(counts))

  if("totalCounts" %in% filters) {
    lqc.totalCounts <- detectLowQualityCells.totalCounts(
      counts,
      mincount = cellMinCount
    )
  }

  if("housekeeping" %in% filters) {
    lqc.housekeeping <- detectLowQualityCells.housekeeping(
      counts,
      geneName = geneName,
      quantileCut = quantileCut.hk
    )
  }

  if("ERCCfrac" %in% filters) {
    lqc.ERCCfrac <- detectLowQualityCells.ERCCfrac(
      counts,
      countsERCC,
      quantileCut = quantileCut.ercc
    )
  }

  lqc <- lqc.totalCounts & lqc.housekeeping & lqc.ERCCfrac

  print(paste0(
    "Removing a total of ",
    sum(!lqc),
    " cells based on the calculated metrics."
  ))

  counts <- counts[, lqc]
  countsERCC <- countsERCC[, lqc]

  #coerce to matrix
  counts <- convertCountsToMatrix(counts)
  countsERCC <- convertCountsToMatrix(countsERCC)

  return(list(counts, countsERCC))
}

#' saveRDA
#'
#' Renames and saves filtered/formated project data as a .rda file in the "data"
#' folder of the current directory.
#'
#' @name saveRDA
#' @rdname saveRDA
#' @param projectName character; The name of the project. Must correspond with
#' the project folder name located at EngeLab/data.
#' @param ... Items to be saved in the .rda file. Note that these variables are
#' renamed before saving such that they have the projectName argument prepended
#' to the current variable name.
#' @return The .rda file is saved as "data/projectName.rda".
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom stringr str_replace

saveRDA <- function(projectName, ...) {
  data <- list(...)
  shortName <- str_replace(projectName, ".*_(.*)", "\\1")
  names <- sapply(substitute(list(...))[-1], deparse)
  names(data) <- paste(shortName, names, sep = ".")

  #reassign the variable name from "counts" etc. to "projectName Counts" etc.
  for(i in 1:length(data)) {
    assign(names(data)[i], data[[i]])
  }

  save(
    list = names(data),
    file = file.path("./data", paste0(projectName, ".rda")),
    compress = "bzip2"
  )
}

#' processedDataUpload
#'
#' Saves the submitted variables as text files in the drives project folder in
#' the "processed_data" subfolder. Note that if the "processed_data" folder is
#' not present in the project dierctory it is created.
#'
#' @name processedDataUpload
#' @rdname processedDataUpload
#' @param projectName character; The name of the project. Must correspond with
#' the project folder name located at EngeLab/data.
#' @param ... Items to be uploaded as .txt files. Note the files are named such
#' that they have the projectName argument prepended to the current variable
#' name.
#' @return The .txt files uploaded to the processed_data folder.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom googledrive drive_ls drive_mkdir drive_upload
#' @importFrom purrr map2 map
#' @importFrom stringr str_replace
#' @importFrom utils write.table

processedDataUpload <- function(projectName, ...) {
  data <- list(...)
  shortName <- str_replace(projectName, ".*_(.*)", "\\1")
  names <- paste(
    shortName,
    sapply(substitute(list(...))[-1], deparse),
    sep = "."
  )
  localPath <- tempdir()

  #check for processed_data dir and create if needed
  subfolders <- googledrive::drive_ls(file.path('data', projectName))$name
  if(!"processed_data" %in% subfolders) {
    googledrive::drive_mkdir("processed_data", file.path('data', projectName))
  }

  #save objects locally
  trash <- purrr::map2(data, names, function(d, n) {
    utils::write.table(
      as.data.frame(d),
      file = file.path(localPath, paste0(n, ".txt")),
      quote = FALSE, sep = "\t"
    )
  })

  #upload to drive
  trash <- purrr::map(names, function(n) {
    googledrive::drive_upload(
      media = file.path(localPath, paste0(n, ".txt")),
      path = file.path(projectName, "processed_data", paste0(n, ".txt"))
    )
  })
  return("Upload successful.")
}
