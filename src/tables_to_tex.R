# read simulation results and transform to latex code
#

library("tidyverse")  # aaaaaah why why why is the default to print all these messages


# 0. preliminaries
# 
# when using this file without the corresponding R-project you have to adjust
# the bld-path.

n_digits <- 3  # number of digits to display in latex tables

r_path <- file.path("bld", "R")  # path to R simulation results
m_path <- file.path("bld", "matlab")  # path to matlab simulation results
tex_path <- file.path("bld", "tex")  # path where to write tex files to


# 1. read data

read_simulations <- function(path) {
  #' Read simulation from path
  #' 
  #' Works for R and matlab simulations
  
  data <- list()
  data$dgp1 <- readr::read_csv(file.path(path, "simulation_dgp1.csv")) %>% mutate_at(c("T", "N"), as.integer)
  tmp_data <- readr::read_csv(file.path(path, "simulation_dgp2-to-dgp4.csv")) %>% mutate_at(c("T", "N", "S"), as.integer) 
  data$dgp2 <- tmp_data %>% filter(dgp == 2)
  data$dgp3 <- tmp_data %>% filter(dgp == 3)
  data$dgp4 <- tmp_data %>% filter(dgp == 4)
  data$dgp5 <- readr::read_csv(file.path(path, "simulation_dgp5.csv")) %>% mutate_at(c("T", "N", "S"), as.integer)
  data$dgp6 <- readr::read_csv(file.path(path, "simulation_dgp6.csv")) %>% mutate_at(c("T", "N"), as.integer)
  
  data[["dgp5"]] <- data[["dgp5"]] %>% add_column(dgp = 5)
  data[["dgp6"]] <- data[["dgp6"]] %>% add_column(dgp = 6, S = 0)
  
  for (dgp in 1:6) {
    data[[paste0("dgp", dgp)]] <- data[[paste0("dgp", dgp)]] %>% round(digits=n_digits)
  }
  
  return(data)
}


data_R <- read_simulations(r_path)
data_m <- read_simulations(m_path)


# 2. create data frames

## 2.2 tables n>30

# remove n = 30 case and irrelevant columns from data

# dgp1 
data_R[["dgp1"]] <- data_R[["dgp1"]] %>% 
  select(-c("s_01", "s_02", "taed_mean", "taed_sd"))
data_m[["dgp1"]] <- data_m[["dgp1"]] %>% 
  select(-c("taed_mean", "taed_sd"))
 
# dgp2 - dgp5
for (dgp in 2:5) {
  data_R[[paste0("dgp", dgp)]] <- data_R[[paste0("dgp", dgp)]] %>% 
    select(-c("dgp", "s_0", "taed_mean", "taed_sd", "mdcj_mean", "mdcj_sd"))
  
  data_m[[paste0("dgp", dgp)]] <- data_m[[paste0("dgp", dgp)]] %>% 
    select(-c("dgp", "s_0", "taed_mean", "taed_sd", "mdcj_mean", "mdcj_sd"))
}

data_R[["dgp5"]] <- data_R[["dgp5"]] %>% select(
  -c("time_effect_mise_mean", "time_effect_mise_sd")
  )


# paper: for the paper we only select all rows corresponding to S = 3 jump loc.

data_R_paper = list()
data_m_paper = list()
for (dgp in 2:5) {
  data_R_paper[[paste0("dgp", dgp)]] <- data_R[[paste0("dgp", dgp)]] %>%
    filter(S == 3)
  
  data_m_paper[[paste0("dgp", dgp)]] <- data_m[[paste0("dgp", dgp)]] %>%
    filter(S == 3)
}

# dgp 6
data_R[["dgp6"]] <- data_R[["dgp6"]] %>% 
  select(-c("dgp", "s_0", "taed_mean", "taed_sd"))

data_m[["dgp6"]] <- data_m[["dgp6"]] %>% 
  select(-c("dgp", "s_0", "taed_mean", "taed_sd"))


# table 1

table1_all <- dplyr::inner_join(
  data_R[["dgp1"]] %>% rename_at(vars(-T, -N), ~ paste0(., ".R")),
  data_m[["dgp1"]] %>% rename_at(vars(-T, -N), ~ paste0(., ".m")),
  by=c("T", "N"),
  suffix=c(".R", ".m")
  )

table1_mise <- table1_all %>%
  select(
    c(
      "T",
      "N",
      "mise1_mean.R",
      "mise1_sd.R",
      "mise2_mean.R",
      "mise2_sd.R",
      "mise1_mean.m",
      "mise1_sd.m",
      "mise2_mean.m",
      "mise2_sd.m"
    )
  )

table1_jumps <- table1_all %>%
  select(
    c(
      "T",
      "N",
      "s_est1_mean.R",
      "s_est1_sd.R",
      "s_est2_mean.R",
      "s_est2_sd.R",
      "hd1_mean.R",
      "hd1_sd.R",
      "hd2_mean.R",
      "hd2_sd.R",
      "s_est_mean.m",
      "s_est_sd.m",
      "hd1_mean.m",
      "hd1_sd.m",
      "hd2_mean.m",
      "hd2_sd.m"
    )
  ) %>% round(digits=2)

table1 = list("all"=table1_all, "mise"=table1_mise, "jumps"=table1_jumps)


# tables 2 - 5

joiner <- function(a, b) {
  out <- dplyr::inner_join(a, b, by=c("T", "S", "N"), suffix=c(".R", ".m"))
  return(out)
}

tables <- lapply(
  paste0("dgp", 2:5),
  function(dgp) joiner(data_R[[dgp]], data_m[[dgp]])
  )

tables_paper <- lapply(
  paste0("dgp", 2:5),
  function(dgp) joiner(data_R_paper[[dgp]], data_m_paper[[dgp]])
  )

tables_paper[[3]] <- tables_paper[[3]] %>% round(digits=2)


# table 6

table6 <- joiner(data_R[["dgp6"]], data_m[["dgp6"]]) %>% select(-c("S"))



################### move sd estimate in bracket behind mean ####################

put_sd_in_brackets <- function(df, index) {
  out <- df %>% select(all_of(index))
  dff = df %>% select(-all_of(index))
  data = lapply(
    seq(1, ncol(dff), 2),
    function(i) do.call(paste0, cbind(dff[i], " (", dff[i+1], ")"))
    )
  cnames = colnames(dff)[seq(1, ncol(dff), 2)]
  data = t(as.data.frame(do.call(rbind, data)))
  colnames(data) = cnames
  out = cbind(out, data)
  rownames(out) = NULL
  return(out) 
}

table1_jumps <- put_sd_in_brackets(table1_jumps, c("T", "N"))
table1_mise <- put_sd_in_brackets(table1_mise, c("T", "N"))

table6 <- put_sd_in_brackets(table6, c("T", "N"))

tables_paper <- lapply(tables_paper, function(df) put_sd_in_brackets(df, c("T", "N", "S")))

tables_supplement <- lapply(tables, function(df) put_sd_in_brackets(df, c("T", "N", "S")))


# 3. Supplement

# For supplement move S to first column and sort rows with respect to (S, T, N)

tables_supplement <- lapply(
  tables_supplement, function(df) df %>% select(S, everything())
  )

tables_supplement <- lapply(
  tables_supplement, function(df) df %>% arrange(S, T, N)
)


############################# tables to latex ##################################

write_table_to_tex <- function(table, fname) {
  table %>%
    xtable::xtable() %>%
    xtable::print.xtable(
      file=file.path(tex_path, paste0(fname, ".tex")), include.rownames = FALSE
      )
}


# table 1

table1_mise %>% write_table_to_tex("table1_mise")
table1_jumps %>% write_table_to_tex("table1_jumps")

# table 2 - 5

for (dgp in 2:5) {
  write_table_to_tex(tables[[dgp - 1]], paste0("table", dgp))
  write_table_to_tex(tables_paper[[dgp - 1]], paste0("table_paper", dgp))
  write_table_to_tex(tables_supplement[[dgp - 1]], paste0("table_supplement", dgp))
}
  
# table 6

table6 %>% write_table_to_tex("table6")

