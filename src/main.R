# main file to run all R simulations

# run SAW monte carlo studies
dir.create(file.path("bld", "R"), showWarnings = FALSE)

source("src/R/monte_carlo_study_dgp1.R")
source("src/R/monte_carlo_study_dgp2.R")
source("src/R/monte_carlo_study_dgp3_to_dgp6.R")

# produce latex tables from simulation results
dir.create(file.path("bld", "tex"), showWarnings = FALSE)
warning("This expects matlab results to be done.")

source("src/tables_to_tex.R")