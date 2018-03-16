#context("metadata_countsMgfp")

library(tibble)
library(dplyr)
data(testingMeta)

test_that("check annotateGFP", {
  
  #setup input data
  names <- c(
    "s.NJB00201.B01", "s.NJB00201.A07", "s.NJB00201.F03", "s.NJB00201.F07",
    "m.NJB00204.A09", "m.NJB00204.B05"
  )
  
  input <- tibble(sample = names) %>%
  annotatePlate(.) %>%
  annotateRow(.) %>%
  annotateColumn(.)
  
  #setup expected data
  expected <- mutate(input, GFP = c(rep(TRUE, 2), rep(FALSE, 2), TRUE, FALSE))
  
  #run function
  output <- annotateGFP(
    input,
    plate = list("NJB00201", "NJB00204"),
    row = list(1:2, 1),
    column = list(c(1, 7), 9)
  )
  
  #test
  expect_identical(expected, output)
})
