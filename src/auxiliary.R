spread_df_over_s <- function(df) {
  tmp_df1 <- df %>% filter(S == 1) %>%
    select(-c(S, T, N)) %>%
    rename_all(tibble::lst( ~ paste0(., '_1')))
  tmp_df2 <- df %>% filter(S == 2) %>%
    select(-c(S, T, N)) %>%
    rename_all(tibble::lst( ~ paste0(., '_2')))
  tmp_df3 <- df %>% filter(S == 3) %>%
    select(-c(S, T, N)) %>%
    rename_all(tibble::lst( ~ paste0(., '_3')))
  
  out <- dplyr::bind_cols(tmp_df1, tmp_df2, tmp_df3)
  return(out)
}

construct_df <- function(df) {
  new_df <- spread_df_over_s(df)
  .id <- df %>% select(c("T", "N")) %>% unique()
  df <- dplyr::bind_cols(.id, new_df)
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


read_simulations <- function(path) {
  #' Read simulation from path
  #' 
  #' Works for R and matlab simulations
  
  data <- list()
  data$dgp1 <- readr::read_csv(file.path(path, "simulation_dgp1.csv"))
  tmp_data <- readr::read_csv(file.path(path, "simulation_dgp2-to-dgp4.csv"))
  data$dgp2 <- tmp_data %>% filter(dgp == 2)
  data$dgp3 <- tmp_data %>% filter(dgp == 3)
  data$dgp4 <- tmp_data %>% filter(dgp == 4)
  data$dgp5 <- readr::read_csv(file.path(path, "simulation_dgp5.csv"))
  data$dgp6 <- readr::read_csv(file.path(path, "simulation_dgp6.csv"))
  
  data[["dgp5"]] <- data[["dgp5"]] %>% add_column(dgp = 5)
  data[["dgp6"]] <- data[["dgp6"]] %>% add_column(dgp = 6, S = 0)
  
  
  return(data)
}

extract_data_for_n30_table <- function(data_list, relevant_columns) {
  
  df <- data_list[["dgp2"]] %>% filter(N == 30) %>% select(all_of(relevant_columns))
  
  for (dgp in 2:6) {
    tmp <- data_list[[paste0("dgp", dgp)]] %>% filter(N == 30) %>% select(all_of(relevant_columns))
    df <- df %>% bind_rows(tmp)
  }
  
  return(df)
}
