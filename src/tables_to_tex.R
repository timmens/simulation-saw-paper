# read simulation results and transform to latex code
#
# the code is structured as follows: first we read in the simulation results
# and perform some basic data cleaning and then we transform the combined
# data frames to latex code

library("magrittr")
library("kableExtra")
library("tidyverse")

# config parameters
n_digits <- 3

r_path <- file.path("bld", "R")
matlab_path <- file.path("bld", "matlab")
tex_path <- file.path("bld", "tex")


############################ auxiliary functions ###############################

spread_df_over_s <- function(df) {
  tmp_df1 <- df %>% filter(S == 1) %>%
    select(-c(S)) %>%
    rename_all(tibble::lst( ~ paste0(., '_1')))
  tmp_df2 <- df %>% filter(S == 2) %>%
    select(-c(S)) %>%
    rename_all(tibble::lst( ~ paste0(., '_2')))
  tmp_df3 <- df %>% filter(S == 3) %>%
    select(-c(S)) %>%
    rename_all(tibble::lst( ~ paste0(., '_3')))
  
  out <- dplyr::bind_cols(tmp_df1, tmp_df2, tmp_df3)
  return(out)
}

construct_df <- function(df) {
  new_df <- spread_df_over_s(df)
  new_cols <- df %>% select(c("T", "N")) %>% distinct()
  
  df <- dplyr::bind_cols(new_cols, new_df)
  return(df)
}

construct_table <- function(table, digits) {
  # does not work for dgp1 and dgp7
  out <-  table %>%
    kable(format = 'latex',
          booktabs = TRUE,
          digits = digits) %>%
    add_header_above(header = c(
      " " = 2,
      "S = 1" = 6,
      "S = 2" = 6,
      "S = 3" = 6
    )) %>%
    kable_styling(latex_options = "striped")
  return(out)
}

########################## tables creation #####################################

## Table 1: data cleaning and construction

# read simulation results
df_R <-
  readr::read_csv(file.path(r_path, "simulation_dgp1.csv")) %>%
  select(-c(s_01, s_02, hd_mean, hd_sd))

df_m <-
  readr::read_csv(file.path(matlab_path, "matlab-simulation-dgp1.csv"))

# add suffix wether results were produces using R or matlab
.id <- df_R[c("T", "N")]

df_R <-
  df_R %>%  mutate(hd1_mean = hd1_mean / T, hd2_mean = hd2_mean / T) %>%
  select(-c(T, N)) %>%
  rename_all(tibble::lst( ~ paste0(., '_R')))

df_m <-
  df_m %>% mutate(hd1_mean = hd1_mean / T, hd2_mean = hd2_mean / T) %>%
  select(-c(T, N)) %>%
  rename_all(tibble::lst( ~ paste0(., '_m')))

# combine R and matlab results
table1 <- dplyr::bind_cols(.id, df_R, df_m)

table1_mise <- table1 %>%
  select(
    c(
      "T",
      "N",
      "mise1_mean_R",
      "mise1_sd_R",
      "mise2_mean_R",
      "mise2_sd_R",
      "mise1_mean_m",
      "mise1_sd_m",
      "mise2_mean_m",
      "mise2_sd_m"
    )
  )

table1_jumps <- table1 %>%
  select(
    c(
      "T",
      "N",
      "s_est1_mean_R",
      "s_est1_sd_R",
      "s_est2_mean_R",
      "s_est2_sd_R",
      "hd1_mean_R",
      "hd1_sd_R",
      "hd2_mean_R",
      "hd2_sd_R",
      "s_est_mean_m",
      "s_est_sd_m",
      "hd1_mean_m",
      "hd1_sd_m",
      "hd2_mean_m",
      "hd2_sd_m"
    )
  )


## Tables 2 - 5 (R): data cleaning and construction for R codes
df <-
  readr::read_csv(file.path(r_path, "simulation_dgp2-to-dgp4.csv")) %>%
  select(-c(mdcj_mean, mdcj_sd, s_0))

df_dgp2 <-
  df %>% filter(dgp == 2) %>% mutate(hd_mean = hd_mean / T) %>% mutate(hd_sd = hd_sd / T)
df_dgp3 <-
  df %>% filter(dgp == 3) %>% mutate(hd_mean = hd_mean / T) %>% mutate(hd_sd = hd_sd / T)
df_dgp4 <-
  df %>% filter(dgp == 4) %>% mutate(hd_mean = hd_mean / T) %>% mutate(hd_sd = hd_sd / T)

table2_R <- construct_df(df_dgp2)
table3_R <- construct_df(df_dgp3)
table4_R <- construct_df(df_dgp4)


## Tables 2 - 4 (matlab): data cleaning and construction for matlab codes

# Update matlab simulationm 

# # # # # # # # # #warning("New simulation go from dgp2-dgp5 not dgp6")
# # # # # # # # # #df <-
# # # # # # # # # #  readr::read_csv(file.path(matlab_path, "matlab-simulation-dgp2-dgp6.csv")) %>%
# # # # # # # # # #  select(-c(mdcj_mean, mdcj_sd, additional_info, s_0))
# # # # # # # # # #
# # # # # # # # # #df_dgp2 <-
# # # # # # # # # #  df %>% filter(dgp == 2) %>% mutate(hd_mean = hd_mean / T) %>% mutate(hd_sd = hd_sd / T)
# # # # # # # # # #df_dgp3 <-
# # # # # # # # # #  df %>% filter(dgp == 3) %>% mutate(hd_mean = hd_mean / T) %>% mutate(hd_sd = hd_sd / T)
# # # # # # # # # #df_dgp4 <-
# # # # # # # # # #  df %>% filter(dgp == 4) %>% mutate(hd_mean = hd_mean / T) %>% mutate(hd_sd = hd_sd / T)
# # # # # # # # # #df_dgp5 <-
# # # # # # # # # #  df %>% filter(dgp == 5) %>% mutate(hd_mean = hd_mean / T) %>% mutate(hd_sd = hd_sd / T)
# # # # # # # # # #
# # # # # # # # # #table2_matlab <- construct_df(df_dgp2)
# # # # # # # # # #table3_matlab <- construct_df(df_dgp3)
# # # # # # # # # #table4_matlab <- construct_df(df_dgp4)
# # # # # # # # # #table5_matlab <- construct_df(df_dgp5)


## Table 5 & 6 (R): data cleaning and construction
table5_R <- readr::read_csv(file.path(r_path, "simulation_dgp5.csv"))
table6_R <- readr::read_csv(file.path(r_path, "simulation_dgp6.csv"))

## Table 6 (matlab): data cleaning and construction

# Update matlba sim

# # # # # df <-
# # # # #   readr::read_csv(file.path(matlabl_path, "matlab-simulation-dgp7.csv"))
# # # # # table7_matlab <- df %>% select(-c("dgp", "additional_info"))



############################# tables to latex ##################################
## R

# create
table1_jumps_tex <- table1_jumps %>%
  kableExtra::kable(
    format = 'latex',
    booktabs = TRUE,
    digits = 2
    ) %>%
  add_header_above(header = c(" " = 2, "R" = 8, "Matlab" = 6)) %>%
  add_header_above(header = c("DGP1" = 16))
# save
kableExtra::save_kable(table1_jumps_tex, file.path(tex_path, "table1_jumps.tex"))


# create
table1_mise_latex <- xtable::xtable(table1_mise, digits = c(0, 0, 0, rep(3, 8)))
# save
xtable::print.xtable(table1_mise_latex, file=file.path(tex_path, "table1_mise.tex"), include.rownames = FALSE)

table2_latex <- construct_table(table2_R, digits = 2) %>%
  add_header_above(header = c("DGP2" = 20))
kableExtra::save_kable(table2_latex, file.path(tex_path, "table2_R.tex"))

table3_latex <- construct_table(table3_R, digits = 2) %>%
  add_header_above(header = c("DGP2" = 20))
kableExtra::save_kable(table3_latex, file.path(tex_path, "table3_R.tex"))

table4_latex <- construct_table(table4_R, digits = 2) %>%
  add_header_above(header = c("DGP2" = 20))
kableExtra::save_kable(table4_latex, file.path(tex_path, "table4_R.tex"))

table5_latex <- construct_table(table5_R, digits=2) %>%
  add_header_above(header = c("DGP2" = 20))
kableExtra::save_kable(table5_latex, file.path(tex_path, "table5_R.tex"))

table6_latex <- table6_R %>%
  kable(format = 'latex',
        booktabs = TRUE,
        digits = n_digits) %>%
  kable_styling(latex_options = "striped") %>%
  add_header_above(header = c("DGP7" = 7))
kableExtra::save_kable(table6_latex, file.path(tex_path, "table6_R.tex"))


## matlab
# # # # # # ## Matlab Simulations
# # # # # # 
# # # # # # construct_table(table2_matlab, digits = n_digits) %>%
# # # # # #   add_header_above(header = c("DGP2" = 20))
# # # # # # 
# # # # # # construct_table(table3_matlab, digits = n_digits) %>%
# # # # # #   add_header_above(header = c("DGP3" = 20))
# # # # # # 
# # # # # # construct_table(table4_matlab, digits = n_digits) %>%
# # # # # #   add_header_above(header = c("DGP4" = 20))
# # # # # # 
# # # # # # construct_table(table5_matlab, digits = n_digits) %>%
# # # # # #   add_header_above(header = c("DGP5" = 20))
# # # # # # 
# # # # # # construct_table(table6_matlab, digits = n_digits) %>%
# # # # # #   add_header_above(header = c("DGP6" = 20))
# # # # # # 
# # # # # # table7_matlab %>%
# # # # # #   kable(format = 'latex',
# # # # # #         booktabs = TRUE,
# # # # # #         digits = n_digits) %>%
# # # # # #   kable_styling(latex_options = "striped") %>%
# # # # # #   add_header_above(header = c("DGP7" = 7))
# # # # # # 
# # # # # # 
# # # # # # # Table 2
# # # # # # # S=1
# # # # # # t2_R_1 <-
# # # # # #   table2 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_1,
# # # # # #                     s_est_sd_1,
# # # # # #                     mise_mean_1,
# # # # # #                     mise_sd_1,
# # # # # #                     hd_mean_1,
# # # # # #                     hd_sd_1)
# # # # # # t2_M_1 <-
# # # # # #   table2_matlab %>% select(s_est_mean_1,
# # # # # #                            s_est_sd_1,
# # # # # #                            mise_mean_1,
# # # # # #                            mise_sd_1,
# # # # # #                            hd_mean_1,
# # # # # #                            hd_sd_1)
# # # # # # t2_1   <- cbind(t2_R_1, t2_M_1)
# # # # # # t2_1_latex <- xtable::xtable(t2_1, digits = c(0, 0, 0, rep(4, 12)))
# # # # # # print(t2_1_latex, include.rownames = FALSE)
# # # # # # # S=2
# # # # # # t2_R_2 <-
# # # # # #   table2 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_2,
# # # # # #                     s_est_sd_2,
# # # # # #                     mise_mean_2,
# # # # # #                     mise_sd_2,
# # # # # #                     hd_mean_2,
# # # # # #                     hd_sd_2)
# # # # # # t2_M_2 <-
# # # # # #   table2_matlab %>% select(s_est_mean_2,
# # # # # #                            s_est_sd_2,
# # # # # #                            mise_mean_2,
# # # # # #                            mise_sd_2,
# # # # # #                            hd_mean_2,
# # # # # #                            hd_sd_2)
# # # # # # t2_2   <- cbind(2, t2_R_2, t2_M_2)
# # # # # # t2_2_latex <- xtable::xtable(t2_2, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t2_2_latex, include.rownames = FALSE)
# # # # # # # S=3
# # # # # # t2_R_3 <-
# # # # # #   table2 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_3,
# # # # # #                     s_est_sd_3,
# # # # # #                     mise_mean_3,
# # # # # #                     mise_sd_3,
# # # # # #                     hd_mean_3,
# # # # # #                     hd_sd_3)
# # # # # # t2_M_3 <-
# # # # # #   table2_matlab %>% select(s_est_mean_3,
# # # # # #                            s_est_sd_3,
# # # # # #                            mise_mean_3,
# # # # # #                            mise_sd_3,
# # # # # #                            hd_mean_3,
# # # # # #                            hd_sd_3)
# # # # # # t2_3   <- cbind(3, t2_R_3, t2_M_3)
# # # # # # t2_3_latex <- xtable::xtable(t2_3, digits = c(0, 0, 0, 0, rep(3, 12)))
# # # # # # print(t2_3_latex, include.rownames = FALSE)
# # # # # # 
# # # # # # 
# # # # # # # Table 3
# # # # # # # S=1
# # # # # # t3_R_1 <-
# # # # # #   table3 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_1,
# # # # # #                     s_est_sd_1,
# # # # # #                     mise_mean_1,
# # # # # #                     mise_sd_1,
# # # # # #                     hd_mean_1,
# # # # # #                     hd_sd_1)
# # # # # # t3_M_1 <-
# # # # # #   table3_matlab %>% select(s_est_mean_1,
# # # # # #                            s_est_sd_1,
# # # # # #                            mise_mean_1,
# # # # # #                            mise_sd_1,
# # # # # #                            hd_mean_1,
# # # # # #                            hd_sd_1)
# # # # # # t3_1   <- cbind(1, t3_R_1, t3_M_1)
# # # # # # t3_1_latex <- xtable::xtable(t3_1, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t3_1_latex, include.rownames = FALSE)
# # # # # # # S=2
# # # # # # t3_R_2 <-
# # # # # #   table3 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_2,
# # # # # #                     s_est_sd_2,
# # # # # #                     mise_mean_2,
# # # # # #                     mise_sd_2,
# # # # # #                     hd_mean_2,
# # # # # #                     hd_sd_2)
# # # # # # t3_M_2 <-
# # # # # #   table3_matlab %>% select(s_est_mean_2,
# # # # # #                            s_est_sd_2,
# # # # # #                            mise_mean_2,
# # # # # #                            mise_sd_2,
# # # # # #                            hd_mean_2,
# # # # # #                            hd_sd_2)
# # # # # # t3_2   <- cbind(2, t3_R_2, t3_M_2)
# # # # # # t3_2_latex <- xtable::xtable(t3_2, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t3_2_latex, include.rownames = FALSE)
# # # # # # # S=3
# # # # # # t3_R_3 <-
# # # # # #   table3 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_3,
# # # # # #                     s_est_sd_3,
# # # # # #                     mise_mean_3,
# # # # # #                     mise_sd_3,
# # # # # #                     hd_mean_3,
# # # # # #                     hd_sd_3)
# # # # # # t3_M_3 <-
# # # # # #   table3_matlab %>% select(s_est_mean_3,
# # # # # #                            s_est_sd_3,
# # # # # #                            mise_mean_3,
# # # # # #                            mise_sd_3,
# # # # # #                            hd_mean_3,
# # # # # #                            hd_sd_3)
# # # # # # t3_3   <- cbind(3, t3_R_3, t3_M_3)
# # # # # # t3_3_latex <- xtable::xtable(t3_3, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t3_3_latex, include.rownames = FALSE)
# # # # # # 
# # # # # # 
# # # # # # # Table 4
# # # # # # # S=1
# # # # # # t4_R_1 <-
# # # # # #   table4 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_1,
# # # # # #                     s_est_sd_1,
# # # # # #                     mise_mean_1,
# # # # # #                     mise_sd_1,
# # # # # #                     hd_mean_1,
# # # # # #                     hd_sd_1)
# # # # # # t4_M_1 <-
# # # # # #   table4_matlab %>% select(s_est_mean_1,
# # # # # #                            s_est_sd_1,
# # # # # #                            mise_mean_1,
# # # # # #                            mise_sd_1,
# # # # # #                            hd_mean_1,
# # # # # #                            hd_sd_1)
# # # # # # t4_1   <- cbind(1, t4_R_1, t4_M_1)
# # # # # # t4_1_latex <- xtable::xtable(t4_1, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t4_1_latex, include.rownames = FALSE)
# # # # # # # S=2
# # # # # # t4_R_2 <-
# # # # # #   table4 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_2,
# # # # # #                     s_est_sd_2,
# # # # # #                     mise_mean_2,
# # # # # #                     mise_sd_2,
# # # # # #                     hd_mean_2,
# # # # # #                     hd_sd_2)
# # # # # # t4_M_2 <-
# # # # # #   table4_matlab %>% select(s_est_mean_2,
# # # # # #                            s_est_sd_2,
# # # # # #                            mise_mean_2,
# # # # # #                            mise_sd_2,
# # # # # #                            hd_mean_2,
# # # # # #                            hd_sd_2)
# # # # # # t4_2   <- cbind(2, t4_R_2, t4_M_2)
# # # # # # t4_2_latex <- xtable::xtable(t4_2, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t4_2_latex, include.rownames = FALSE)
# # # # # # # S=3
# # # # # # t4_R_3 <-
# # # # # #   table4 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_3,
# # # # # #                     s_est_sd_3,
# # # # # #                     mise_mean_3,
# # # # # #                     mise_sd_3,
# # # # # #                     hd_mean_3,
# # # # # #                     hd_sd_3)
# # # # # # t4_M_3 <-
# # # # # #   table4_matlab %>% select(s_est_mean_3,
# # # # # #                            s_est_sd_3,
# # # # # #                            mise_mean_3,
# # # # # #                            mise_sd_3,
# # # # # #                            hd_mean_3,
# # # # # #                            hd_sd_3)
# # # # # # t4_3   <- cbind(3, t4_R_3, t4_M_3)
# # # # # # t4_3_latex <- xtable::xtable(t4_3, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t4_3_latex, include.rownames = FALSE)
# # # # # # 
# # # # # # 
# # # # # # # Table 5
# # # # # # # S=1
# # # # # # t5_R_1 <-
# # # # # #   table5 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_1,
# # # # # #                     s_est_sd_1,
# # # # # #                     mise_mean_1,
# # # # # #                     mise_sd_1,
# # # # # #                     hd_mean_1,
# # # # # #                     hd_sd_1)
# # # # # # t5_M_1 <-
# # # # # #   table5_matlab %>% select(s_est_mean_1,
# # # # # #                            s_est_sd_1,
# # # # # #                            mise_mean_1,
# # # # # #                            mise_sd_1,
# # # # # #                            hd_mean_1,
# # # # # #                            hd_sd_1)
# # # # # # t5_1   <- cbind(1, t5_R_1, t5_M_1)
# # # # # # t5_1_latex <- xtable::xtable(t5_1, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t5_1_latex, include.rownames = FALSE)
# # # # # # # S=2
# # # # # # t5_R_2 <-
# # # # # #   table5 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_2,
# # # # # #                     s_est_sd_2,
# # # # # #                     mise_mean_2,
# # # # # #                     mise_sd_2,
# # # # # #                     hd_mean_2,
# # # # # #                     hd_sd_2)
# # # # # # t5_M_2 <-
# # # # # #   table5_matlab %>% select(s_est_mean_2,
# # # # # #                            s_est_sd_2,
# # # # # #                            mise_mean_2,
# # # # # #                            mise_sd_2,
# # # # # #                            hd_mean_2,
# # # # # #                            hd_sd_2)
# # # # # # t5_2   <- cbind(2, t5_R_2, t5_M_2)
# # # # # # t5_2_latex <- xtable::xtable(t5_2, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t5_2_latex, include.rownames = FALSE)
# # # # # # # S=3
# # # # # # t5_R_3 <-
# # # # # #   table5 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_3,
# # # # # #                     s_est_sd_3,
# # # # # #                     mise_mean_3,
# # # # # #                     mise_sd_3,
# # # # # #                     hd_mean_3,
# # # # # #                     hd_sd_3)
# # # # # # t5_M_3 <-
# # # # # #   table5_matlab %>% select(s_est_mean_3,
# # # # # #                            s_est_sd_3,
# # # # # #                            mise_mean_3,
# # # # # #                            mise_sd_3,
# # # # # #                            hd_mean_3,
# # # # # #                            hd_sd_3)
# # # # # # t5_3   <- cbind(3, t5_R_3, t5_M_3)
# # # # # # t5_3_latex <- xtable::xtable(t5_3, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t5_3_latex, include.rownames = FALSE)
# # # # # # 
# # # # # # 
# # # # # # 
# # # # # # # Table 6
# # # # # # # S=1
# # # # # # t6_R_1 <-
# # # # # #   table6 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_1,
# # # # # #                     s_est_sd_1,
# # # # # #                     mise_mean_1,
# # # # # #                     mise_sd_1,
# # # # # #                     hd_mean_1,
# # # # # #                     hd_sd_1)
# # # # # # t6_M_1 <-
# # # # # #   table6_matlab %>% select(s_est_mean_1,
# # # # # #                            s_est_sd_1,
# # # # # #                            mise_mean_1,
# # # # # #                            mise_sd_1,
# # # # # #                            hd_mean_1,
# # # # # #                            hd_sd_1)
# # # # # # t6_1   <- cbind(1, t6_R_1, t6_M_1)
# # # # # # t6_1_latex <- xtable::xtable(t6_1, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t6_1_latex, include.rownames = FALSE)
# # # # # # # S=2
# # # # # # t6_R_2 <-
# # # # # #   table6 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_2,
# # # # # #                     s_est_sd_2,
# # # # # #                     mise_mean_2,
# # # # # #                     mise_sd_2,
# # # # # #                     hd_mean_2,
# # # # # #                     hd_sd_2)
# # # # # # t6_M_2 <-
# # # # # #   table6_matlab %>% select(s_est_mean_2,
# # # # # #                            s_est_sd_2,
# # # # # #                            mise_mean_2,
# # # # # #                            mise_sd_2,
# # # # # #                            hd_mean_2,
# # # # # #                            hd_sd_2)
# # # # # # t6_2   <- cbind(2, t6_R_2, t6_M_2)
# # # # # # t6_2_latex <- xtable::xtable(t6_2, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t6_2_latex, include.rownames = FALSE)
# # # # # # # S=3
# # # # # # t6_R_3 <-
# # # # # #   table6 %>% select(T,
# # # # # #                     N,
# # # # # #                     s_est_mean_3,
# # # # # #                     s_est_sd_3,
# # # # # #                     mise_mean_3,
# # # # # #                     mise_sd_3,
# # # # # #                     hd_mean_3,
# # # # # #                     hd_sd_3)
# # # # # # t6_M_3 <-
# # # # # #   table6_matlab %>% select(s_est_mean_3,
# # # # # #                            s_est_sd_3,
# # # # # #                            mise_mean_3,
# # # # # #                            mise_sd_3,
# # # # # #                            hd_mean_3,
# # # # # #                            hd_sd_3)
# # # # # # t6_3   <- cbind(3, t6_R_3, t6_M_3)
# # # # # # t6_3_latex <- xtable::xtable(t6_3, digits = c(0, 0, 0, 0, rep(2, 12)))
# # # # # # print(t6_3_latex, include.rownames = FALSE)
# # # # # # 
# # # # # # 
# # # # # # 
# # # # # # # Table 7
# # # # # # # S=0
# # # # # # t7_R <-
# # # # # #   table7 %>% select(T, N,  s_est_mean, s_est_sd, mise_mean, mise_sd)
# # # # # # t7_M <-
# # # # # #   table7_matlab %>% select(s_est_mean, s_est_sd, mise_mean, mise_sd)
# # # # # # t7   <- cbind(0, t7_R, t7_M)
# # # # # # t7_latex <- xtable::xtable(t7, digits = c(0, 0, 0, 0, 2, 2, 3, 3, 2, 2, 3, 3))
# # # # # # print(t7_latex, include.rownames = FALSE)