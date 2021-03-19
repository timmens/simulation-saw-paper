# main file to run all R simulations

output_dir = file.path("bld", "R")
dir.create(file.path(output_dir), showWarnings = FALSE)

# run monte carlo studies
source("src/R/monte_carlo_study_dgp1.R")
source("src/R/monte_carlo_study_dgp2.R")
source("src/R/monte_carlo_study_dgp3_to_dgp6.R")

# produce extra output given simulation results