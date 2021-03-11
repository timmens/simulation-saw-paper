library(magrittr)
library(kableExtra)
library(tidyverse)

spread_df_over_s <- function(df) {
  tmp_df <- df %>% select(-c(1, 2, 4))
  
  tmp_df1 <- tmp_df %>% filter(S == 1) %>% 
    select(-c(S)) %>%
    rename_all(tibble::lst(~ paste0(., '_1')))
  tmp_df2 <- tmp_df %>% filter(S == 2) %>% 
    select(-c(S)) %>%
    rename_all(tibble::lst(~ paste0(., '_2')))
  tmp_df3 <- tmp_df %>% filter(S == 3) %>% 
    select(-c(S)) %>%
    rename_all(tibble::lst(~ paste0(., '_3')))
  
  dplyr::bind_cols(tmp_df1, tmp_df2, tmp_df3)
}

construct_df <- function(df) {
  new_df <- spread_df_over_s(df)
  new_cols <- df %>% select(c("T", "N")) %>% distinct()
  df <- dplyr::bind_cols(new_cols, new_df)
}

construct_table <- function(table, digits) {
  # does not work for dgp1 and dgp7 
  out <-  table %>%
    kable(format = 'latex', booktabs=TRUE, digits=digits) %>%
    add_header_above(header = c(" " = 2, "S = 1" = 3, "S = 2" = 3, "S = 3" = 3)) %>%
    kable_styling(latex_options = "striped")
  out
}

### Table for DGP 1
df_R <- readr::read_csv("rsimulation-dgp1-alt-2019-11-04-17-05-28.csv")
df_m <- readr::read_csv("../../matlab_dgp1_simulation.csv") %>% 
  select(-c(additional_info)) %>% rename(mise1 = mse1) %>% rename(mise2 = mse2)

df_ind <- df_R[,1:2]
df_R <- df_R %>%  mutate(hd1 = hd1 / T, hd2 = hd2 / T) %>% select(-c(T, N)) %>% 
  rename_all(tibble::lst(~ paste0(., '_R')))
df_m <- df_m %>% mutate(hd1 = hd1 / T, hd2 = hd2 / T) %>% select(-c(T, N)) %>%
  rename_all(tibble::lst(~ paste0(., '_m')))

table1 <- dplyr::bind_cols(df_ind, df_R, df_m)

table1 %>% 
  kable(format = 'latex', booktabs=TRUE, digits=3) %>%
  kable_styling(latex_options = "striped") %>% 
  add_header_above(header = c(" " = 2, "R" = 4, "Matlab" = 3)) %>% 
  add_header_above(header = c("DGP1" = 9))
