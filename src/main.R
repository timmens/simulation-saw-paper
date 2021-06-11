# main file to run all R simulations

# install sawr packagae if necessary
# library("devtools")
# devtools::install_github("https://github.com/timmens/sawr", force = TRUE)

dir.create(file.path("bld", "R"), showWarnings = FALSE)
dir.create(file.path("bld", "tex"), showWarnings = FALSE)

run_simulations <- TRUE

if (run_simulations) {
  
  source("src/R/monte_carlo_study_dgp1.R")
  source("src/R/monte_carlo_study_dgp2_to_dgp4.R")
  source("src/R/monte_carlo_study_dgp5.R")
  source("src/R/monte_carlo_study_dgp6.R")
  
}

source("src/tables_to_tex.R")
