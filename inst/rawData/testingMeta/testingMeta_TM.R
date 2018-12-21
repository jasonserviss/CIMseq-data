#test data for metadata functions
#run with: source('./inst/rawData/testingMeta/testingMeta.R')

TM <- function(upload = FALSE, save = TRUE) {
  cat('Processing testingMeta.\n')

  testingMeta <- tibble::tibble(
    sample = c(
      "s.NJB00201.B01", "s.NJB00201.A07", "s.NJB00201.F03", "s.NJB00201.F07",
      "m.NJB00204.A09", "m.NJB00204.B05")
  )

  if(save) save(
    testingMeta, file = "./data/testingMeta.rda", compress = "bzip2"
  )
}
