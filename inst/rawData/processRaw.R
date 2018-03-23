
ignore <- c(
  "./inst/rawData//permutations_171116/permutations.R",
  "./inst/rawData//syntheticData/syntheticData.R"
)

toProcess <- list.files(
  path = './inst/rawData/', recursive = TRUE,
  pattern = '\\.R$', full.names = TRUE
)

toProcess <- toProcess[!toProcess %in% ignore]

purrr::map(toProcess, source)
