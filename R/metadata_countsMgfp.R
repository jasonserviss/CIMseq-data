
#' Load meta data.
#'
#' Loads a metadata excel file into R.
#'
#' @name loadMetaData
#' @rdname loadMetaData
#' @param path character; The path to the metadata file. Typically in a
#'  subfolder of inst/rawData.
#' @return The metadata as a tibble.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom readxl read_excel

loadMetaData <- function(path) {
  read_excel(path)
}

#' Annotate plate.
#'
#' Uses the standard sample naming nomenclature to add plate row to metadata.
#' Sample names should follow: (s|m)\\.platename\\.platePosition.
#' platePosition is in the form row and column without a space where row is a
#' LETTER (A-H) and column is a number (1-12).
#'
#' @name annotatePlate
#' @rdname annotatePlate
#' @param data tibble; A tibble containing the sample names using standard
#'  nomenclature in a column named "sample".
#' @return The metadata tibble with an additional column indicating plate name.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace

annotatePlate <- function(data) {
  mutate(data, plate = str_replace(sample, "^.\\.([A-Z0-9]*)\\....", "\\1"))
}

#' Annotate row.
#'
#' Uses the standard sample naming nomenclature to add plate row to metadata.
#' Sample names should follow: (s|m)\\.platename\\.platePosition.
#' platePosition is in the form row and column without a space where row is a
#' LETTER (A-H) and column is a number (1-12).
#'
#' @name annotateRow
#' @rdname annotateRow
#' @param data tibble; A tibble containing the sample names using standard
#'  nomenclature in a column named "sample".
#' @return The metadata tibble with an additional column indicating row as a
#' numeric value.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate select
#' @importFrom stringr str_replace

annotateRow <- function(data) {
  data %>%
  mutate(rowPos = str_replace(sample, "^.\\.[A-Z0-9]*\\.(.)..", "\\1")) %>%
  mutate(row = match(rowPos, LETTERS[1:8])) %>%
  select(-rowPos)
}

#' Annotate column.
#'
#' Uses the standard sample naming nomenclature to add plate row to metadata.
#' Sample names should follow: (s|m)\\.platename\\.platePosition.
#' platePosition is in the form row and column without a space where row is a
#' LETTER (A-H) and column is a number (1-12).
#'
#' @name annotateColumn
#' @rdname annotateColumn
#' @param data tibble; A tibble containing the sample names using standard
#'  nomenclature in a column named "sample".
#' @return The metadata tibble with an additional column indicating column as a
#' numeric value.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate case_when select
#' @importFrom stringr str_replace

annotateColumn <- function(data) {
  data %>%
  mutate(colPos = str_replace(sample, "^.\\.[A-Z0-9]*\\..(..)", "\\1")) %>%
  mutate(column = case_when(
   colPos == "01" ~ 1L, colPos == "02" ~ 2L, colPos == "03" ~ 3L,
   colPos == "04" ~ 4L, colPos == "05" ~ 5L, colPos == "06" ~ 6L,
   colPos == "07" ~ 7L, colPos == "08" ~ 8L, colPos == "09" ~ 9L,
   colPos == "10" ~ 10L, colPos == "11" ~ 11L, colPos == "12" ~ 12L,
   colPos == "13" ~ 13L, colPos == "14" ~ 14L, colPos == "15" ~ 15L,
   colPos == "16" ~ 16L, colPos == "17" ~ 17L, colPos == "18" ~ 18L
  )) %>%
  select(-colPos)
}

#' Annotate sample GFP status.
#'
#' Annotates a logical value indicating a sample's GFP status.
#'
#' @name annotateGFP
#' @rdname annotateGFP
#' @param data tibble; A tibble containing the plate, row, and column
#'  designation for the samples to be annotated.
#' @param plate character; A character string indicating the plate variable
#'  where samples to be annotated are located.
#' @param row numeric; A numeric vector indicating the rows that are GFP
#'  positive.
#' @param column numeric; A numeric vector indicating the columns that are GFP
#'  positive.
#' @return The metadata tibble with an additional column indicating GFP staus.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate if_else

annotateGFP <- function(data, plate, row, column) {
  #rename args
  plateIn <- plate; rowIn <- row; colIn <- column
  
  #process
  data %>%
  mutate(
    GFP = if_else(
      plate %in% plateIn & row %in% rowIn & column %in% colIn,
      TRUE,
      FALSE
    )
  )
}
