# read simulation results and transform to latex code
#

library("magrittr")
library("kableExtra")
library("tidyverse", quietly = TRUE)  # aaaaaah why why why is the default to print all these messages

source("src/auxiliary.R")
# exports: `read_simulations`, `extract_data_for_n30_table`, `construct_df`,
# `construct_table`


# 0. preliminaries

n_digits <- 3  # number of digits to display in latex tables

r_path <- file.path("bld", "R")
m_path <- file.path("bld", "matlab")
tex_path <- file.path("bld", "tex")


# 1. read data

data_R <- read_simulations(r_path)
data_m <- read_simulations(m_path)


# 2. create data frames


## 2.1 table n=30

n30_columns <- c("dgp", "T", "S", "s_est_mean", "s_est_sd", "taed_mean", "taed_sd")

df_n30_R <- extract_data_for_n30_table(data_R, n30_columns)
df_n30_m <- extract_data_for_n30_table(data_m, n30_columns)

table_n30 <- dplyr::inner_join(df_n30_R, df_n30_m, by=c("dgp", "T", "S"), suffix=c(".R", ".m"))

## 2.2 tables n>30

# remove n = 30 case and irrelevant columns from data

# dgp1 
data_R[["dgp1"]] <- data_R[["dgp1"]] %>% 
  filter(N != 30) %>% 
  select(-c("s_01", "s_02", "taed_mean", "taed_sd"))
data_m[["dgp1"]] <- data_m[["dgp1"]] %>% 
  filter(N != 30) %>% 
  select(-c("taed_mean", "taed_sd"))
 
# dgp2 - dgp5
for (dgp in 2:5) {
  data_R[[paste0("dgp", dgp)]] <- data_R[[paste0("dgp", dgp)]] %>% 
    filter(N != 30) %>% 
    select(-c("dgp", "s_0", "taed_mean", "taed_sd", "mdcj_mean", "mdcj_sd"))
  
  data_m[[paste0("dgp", dgp)]] <- data_m[[paste0("dgp", dgp)]] %>% 
    filter(N != 30) %>% 
    select(-c("dgp", "s_0", "taed_mean", "taed_sd", "mdcj_mean", "mdcj_sd"))
}

# dgp 6
data_R[["dgp6"]] <- data_R[["dgp6"]] %>% 
  filter(N != 30) %>% 
  select(-c("dgp", "s_0", "taed_mean", "taed_sd"))

data_m[["dgp6"]] <- data_m[["dgp6"]] %>% 
  filter(N != 30) %>% 
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
  )

table1 = list("all"=table1_all, "mise"=table1_mise, "jumps"=table1_jumps)


# tables 2 - 5

table2_R <- construct_df(data_R[["dgp2"]])
table3_R <- construct_df(data_R[["dgp3"]])
table4_R <- construct_df(data_R[["dgp4"]])
table5_R <- construct_df(
  data_R[["dgp5"]] %>% select(-c("time_effect_mise_mean", "time_effect_mise_sd"))
  )

table2_m <- construct_df(data_m[["dgp2"]])
table3_m <- construct_df(data_m[["dgp3"]])
table4_m <- construct_df(data_m[["dgp4"]])
table5_m <- construct_df(data_m[["dgp5"]])


# table 5: time fixed effect part

table5_time_effect <- data_R[["dgp5"]] %>% select(
  -c(s_est_mean, s_est_sd, mise_mean, mise_sd, hd_mean, hd_sd)
  )

warning("Look at the data one more time, sd needs to be non-zero")


# table 6

table6 <- dplyr::inner_join(
  data_R[["dgp6"]],
  data_m[["dgp6"]], 
  by=c("T", "N", "S"),
  suffix=c(".R", ".m")
  ) %>% select(-c("S"))



############################# tables to latex ##################################

# table 1

## create

table1_jumps_tex <- table1_jumps %>%
  kableExtra::kable(
    format = 'latex',
    booktabs = TRUE,
    digits = 2
    ) %>%
  add_header_above(header = c(" " = 2, "R" = 8, "Matlab" = 6)) %>%
  add_header_above(header = c("DGP1" = 16))

table1_mise_latex <- xtable::xtable(table1_mise, digits = c(0, 0, 0, rep(3, 8)))

## save

kableExtra::save_kable(table1_jumps_tex, file.path(tex_path, "table1_jumps.tex"))

xtable::print.xtable(table1_mise_latex, file=file.path(tex_path, "table1_mise.tex"), include.rownames = FALSE)


# table 6

## create

table6_latex <- table6 %>%
  kable(format = 'latex',
        booktabs = TRUE,
        digits = n_digits) %>%
  kable_styling(latex_options = "striped") %>%
  add_header_above(header = c("DGP7" = 7))

## save

kableExtra::save_kable(table6_latex, file.path(tex_path, "table6_R.tex"))



##### UNCLEAR WHAT HAPPENED BELOW HERE

# # table2_latex <- construct_table(table2_R, digits = 2) %>%
# #   add_header_above(header = c("DGP2" = 20))
# # kableExtra::save_kable(table2_latex, file.path(tex_path, "table2_R.tex"))
# # 
# # table3_latex <- construct_table(table3_R, digits = 2) %>%
# #   add_header_above(header = c("DGP2" = 20))
# # kableExtra::save_kable(table3_latex, file.path(tex_path, "table3_R.tex"))
# # 
# # table4_latex <- construct_table(table4_R, digits = 2) %>%
# #   add_header_above(header = c("DGP2" = 20))
# # kableExtra::save_kable(table4_latex, file.path(tex_path, "table4_R.tex"))
# # 
# # table5_latex <- construct_table(table5_R, digits=2) %>%
# #   add_header_above(header = c("DGP2" = 20))
# # kableExtra::save_kable(table5_latex, file.path(tex_path, "table5_R.tex"))



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