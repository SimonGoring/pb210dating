library(dplyr)
library(pb210dating)

raw_data <- readr::read_csv('inst/ABRA2 20171006.csv')

pb_const <- constants(file = "filename.txt", raw_data)

# Estimate visually:
plot(pb_const)

pb_const <- constants(file = "filename.txt", raw_data,
                      "Pb210", "Ra226", NSectionDating = 12)

cfcs()
