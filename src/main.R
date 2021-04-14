# main file to run all R simulations

dir.create(file.path("bld", "R"), showWarnings = FALSE)
dir.create(file.path("bld", "tex"), showWarnings = FALSE)

# run SAW monte carlo studies

run_simulations <- readline("Should simulations be run again?\nType 'yes' for yes and press ENTER for no: ")
run_simulations <- ifelse(run_simulations == "yes", TRUE, FALSE)

if (run_simulations) {
  
  
  source("src/R/monte_carlo_study_dgp1.R")
  source("src/R/monte_carlo_study_dgp2_to_dgp4.R")
  source("src/R/monte_carlo_study_dgp5.R")
  source("src/R/monte_carlo_study_dgp6.R")
  
}

source("src/tables_to_tex.R")
