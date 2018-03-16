
#' Annotate cell types in each well for the countsSorted2 data.
#'
#'
#' @name countsSorted2_cellTypes
#' @rdname countsSorted2_cellTypes
#' @param data tibble; A tibble containing the sample names using standard
#'  nomenclature in a column named "sample".
#' @return The metadata tibble with an additional column, named "cellTypes",
#'  indicating the cell types in the sample. In the case of multiplets, cell
#'  types are seperated by "-".
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate if_else case_when full_join arrange pull n
#' @importFrom tibble as_tibble
#' @importFrom stringr str_extract

countsSorted2_cellTypes <- function(data) {
  data %>%
  mutate(cellTypes = if_else(
    cellNumber == "Singlet",
    .annotateCellTypeSng(sample),
    .annotateCellTypeMul(sample)
  ))

}

.annotateCellTypeSng <- function(sample) {
  positions <- str_extract(sample, "...$")
  case_when(
    positions == "E03" ~ "A375", #what is this? sorting issue?
    positions %in% paste0(sort(rep(LETTERS[1:8], 4)), c("01", "02", "03", "04")) ~ "HOS",
    positions %in% paste0(sort(rep(LETTERS[1:8], 4)), c("05", "06", "07", "08")) ~ "HCT116",
    positions %in% paste0(sort(rep(LETTERS[1:8], 4)), c("09", "10", "11", "12")) ~ "A375",
    TRUE ~ "error"
  )
}

.annotateCellTypeMul <- function(sample) {
  
  #make plate data
  multiplets <- c(
    c(rep("A375-HCT116", 3), rep("HCT116-HOS", 2), rep("A375-HOS", 3)),
    c(rep("A375-HCT116", 3), rep("HCT116-HOS", 2), rep("A375-HOS", 3)),
    c(rep("A375-HCT116", 2), rep("HCT116-HOS", 3), rep("A375-HOS", 3)),
    c(rep("A375-HCT116", 2), rep("HCT116-HOS", 3), rep("A375-HOS", 3)),
    c(rep("A375-HCT116-HOS", 3), rep("HCT116-HCT116-HOS", 2), rep("A375-HOS-HOS", 3)),
    c(rep("A375-HCT116-HOS", 3), rep("HCT116-HCT116-HOS", 2), rep("A375-HOS-HOS", 3)),
    c(rep("A375-HCT116-HOS", 2), rep("HCT116-HCT116-HOS", 3), rep("A375-HOS-HOS", 3)),
    c(rep("A375-HCT116-HOS", 2), rep("HCT116-HCT116-HOS", 3), rep("A375-HOS-HOS", 3)),
    c(rep("A375-HOS-HOS-HOS", 3), rep("A375-A375-HCT116-HCT116", 2), rep("A375-A375-HCT116-HOS", 3)),
    c(rep("A375-HOS-HOS-HOS", 3), rep("A375-A375-HCT116-HCT116", 2), rep("A375-A375-HCT116-HOS", 3)),
    c(rep("A375-HOS-HOS-HOS", 2), rep("A375-A375-HCT116-HCT116", 3), rep("A375-A375-HCT116-HOS", 3)),
    c(rep("A375-HOS-HOS-HOS", 2), rep("A375-A375-HCT116-HCT116", 3), rep("A375-A375-HCT116-HOS", 3))
  )
  
  cols <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  rows <- LETTERS[1:8]
  names <- paste0("m.NJB00204.", rep(rows, 12), sort(rep(cols, 8)))
  
  plateData <- tibble(
    value = names,
    cellTypes = multiplets
  )
  
  sample %>%
  as_tibble() %>%
  mutate(order = 1:n()) %>%
  full_join(plateData, by = "value") %>%
  arrange(order) %>%
  pull(cellTypes)
}







