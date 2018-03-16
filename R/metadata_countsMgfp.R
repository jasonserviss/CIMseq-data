
#' Annotate sample GFP status.
#'
#' Annotates a logical value indicating a sample's GFP status.
#'
#' @name annotateGFP
#' @rdname annotateGFP
#' @param data tibble; A tibble containing the plate, row, and column
#'  designation for the samples to be annotated.
#' @param plate list; A list where each element is a character vector indicating
#'  the plate where samples to be annotated are located.
#' @param row list; A list where each element is a numeric vector indicating the
#'  rows that are GFP positive in the corresponding plate.
#' @param column list; A list where each element is a numeric vector indicating
#'  the columns that are GFP positive in the corresponding plate.
#' @return The metadata tibble with an additional column indicating GFP staus.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate if_else
#' @importFrom purrr pmap

annotateGFP <- function(data, plate, row, column) {
  
  bool <- pmap(list(plate, row, column), function(x, y, z) {
    data$plate %in% x & data$row %in% y & data$column %in% z
  }) %>%
  Reduce("+", .) %>%
  as.logical()
  
  mutate(data, GFP = if_else(bool, TRUE, FALSE))
}

#' Annotate mouse.
#'
#' Annotates the mouse ID for the experiment. For the moment, assumes that each
#' plate contains only one mouse.
#'
#' @name annotateMouse
#' @rdname annotateMouse
#' @param data tibble; A tibble containing the plate designation for the samples
#'  to be annotated.
#' @param plate character; A character vector indicating the plate where samples
#'  to be annotated are located.
#' @param mouseID character; ID of the mouse for the corresponding plate.
#' @return The metadata tibble with an additional column indicating mouse ID in
#'  a column named "mouseID".
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate if_else

annotateMouse <- function(data, plate, mouseID) {
  plateIN <- plate; mouseID.IN <- mouseID
  mutate(data, mouseID = mouseID.IN[match(plate, plateIN)])
}

#' Annotate tissue.
#'
#' Annotates the tissue for the experiment. For the moment, assumes that each
#' plate contains only one tissue.
#'
#' @name annotateTissue
#' @rdname annotateTissue
#' @param data tibble; A tibble containing the plate designation for the samples
#'  to be annotated.
#' @param plate character; A character vector indicating the plate where samples
#'  to be annotated are located.
#' @param tissue character; Tissue for the corresponding plate.
#' @return The metadata tibble with an additional column indicating tissue in
#'  a column named "tissue".
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate if_else

annotateTissue <- function(data, plate, tissue) {
  plateIN <- plate; tissue.IN <- tissue
  mutate(data, tissue = tissue.IN[match(plate, plateIN)])
}

#' Rename GFP mouse samples.
#'
#' Updates sample names from the old nomenclature (GFPpos.C1.Doublet.5.F9.htseq)
#' to the new nomenclature (NJA00103.D09.htseq).
#'
#' @name renameMgfpSamples
#' @rdname renameMgfpSamples
#' @param oldNames character; A character vector of the old names that have
#'  already been processed by the labelSingletsAndMultiplets and
#'  removeHTSEQsuffix functions. Names that are already in the new nomenclature
#'  will not be effected.
#' @return A character vector of the updated names.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom stringr str_detect str_replace str_extract
#' @importFrom dplyr case_when

renameMgfpSamples <- function(oldNames) {
  bool1 <- str_detect(oldNames, "GFPneg.C1.Singlet.4")
  bool2 <- str_detect(oldNames, "GFPneg.SI1.Singlet.2")
  bool3 <- str_detect(oldNames, "GFPpos.C1.Doublet.5")
  bool4 <- str_detect(oldNames, "GFPpos.C1.Singlet.4")
  bool5 <- str_detect(oldNames, "GFPpos.SI1.Doublet.3")
  bool <- bool1 | bool2 | bool3 | bool4 | bool5
  idx <- which(bool)
  
  #extract and fix plate positions
  pos <- paste0(".", str_extract(oldNames[bool], "[^.]+$"))
  posidx <- which(nchar(pos) == 3)
  pos[posidx] <- str_replace(pos[posidx], "(..)(.)", "\\10\\2")
  
  #extract prefix
  prefix <- str_extract(oldNames[bool], "^..")
  
  #update plate info
  plate <- case_when(
    bool1 ~ "NJA00110",
    bool2 ~ "NJA00101",
    bool3 ~ "NJA00111",
    bool4 ~ "NJA00110",
    bool5 ~ "NJA00107",
    TRUE ~ oldNames
  )
  
  #synthesize new names
  plate[bool] <- paste0(prefix, plate[bool], pos)
  plate
}
