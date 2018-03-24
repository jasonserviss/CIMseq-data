#source('./inst/rawData/processRaw.R')

#run scripts to generate data/.rda files
ignore <- c(
  './inst/rawData/permutations_171116/permutations.R',
  './inst/rawData/syntheticData/syntheticData.R',
  './inst/rawData/processRaw.R'
)

toProcess <- list.files(
  path = './inst/rawData', recursive = TRUE,
  pattern = '\\.R$', full.names = TRUE
)

toProcess <- toProcess[!toProcess %in% ignore]

purrr::map(toProcess, source)

#check output
generateChecksum <- function(filePath) {
  tools::md5sum(filePath)
}

expected <- c(
  "./data/countsMgfp.rda" = "e22884946e164a66f9fa19f67bbddceb",
  "./data/countsRegev.rda" = "74eca56fddc5f041e2dc96ba771ca2a1",
  "./data/countsSorted1.rda" = "d406e4cb5d2e92095b407cd72f5bbe9d",
  "./data/countsSorted2.rda" = "e7f799f0481409073fd58da632141c00",
  "./data/fetalPancreasCounts.rda" = "b2875460d6ac5fbcf116c32308ad5cb0",
  "./data/testingCounts.rda" = "4af88bfc027781f40a4b4e6e01f7220f",
  "./data/testingMeta.rda" = "95f895f9a6f7eff0a85c867af816a1c4"
)
detected <- tools::md5sum(names(expected))
if(!identical(detected, expected)) {
  problem <- names(expected)[which(expected != detected)]
  message <- paste0("Checksums did not match for file(s) ", problem)
}
