#' Rename GFP mouse samples.
#'
#' Updates sample names from the old nomenclature (GFPpos.C1.Doublet.5.F9.htseq)
#' to the new nomenclature (NJA00103.D09.htseq).
#'
#' @name renameMgfpSamples
#' @rdname renameMgfpSamples
#' @param oldNames character; A character vector of the old names that have
#'  already been processed by the removeHTSEQsuffix functions. Names that are
#' already in the new nomenclature will not be effected.
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
  plate[bool] <- paste0(plate[bool], pos)
  plate
}
