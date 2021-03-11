---
title: LaTeX Tables
author: Tim Mensinger
output: pdf_document
--- 

```{r}
library(magrittr)
library(kableExtra)
library(tidyverse)

### Tables for DGP 2 - DGP 4 

df <- readr::read_csv("FINAL_SIM_RESULTS/rsimulaton-dgp2-to-dgp4-final.csv") %>% select(-c(1))
df <- df %>% select(-c(mdcj, additional_info))

df_dgp2 <- df %>% filter(DGP == 2)
df_dgp3 <- df %>% filter(DGP == 3)
df_dgp4 <- df %>% filter(DGP == 4)

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

table2 <- construct_df(df_dgp2)
table3 <- construct_df(df_dgp3)
table4 <- construct_df(df_dgp4)

### Tables for DGP 5 - DGP 6

fd <- readr::read_csv("FINAL_SIM_RESULTS/rsimulaton-dgp5-to-dgp6-final.csv") %>% select(-c(1))
fd <- fd %>% select(-c(mdcj, additional_info))

df_dgp5 <- fd %>% filter(DGP == 5)
df_dgp6 <- fd %>% filter(DGP == 6)

table5 <- construct_df(df_dgp5)
table6 <- construct_df(df_dgp6)
```

```{r eval=FALSE}
### Construct LaTeX Tables 

table2 %>%
    kable(format = 'latex', booktabs = TRUE) %>%
    add_header_above(header = c("Text" = 2, "Values" = 2))
```



\begin{tabular}{rrrrrrrrrrrrrr}
\toprule
\multicolumn{2}{c}{Text} & \multicolumn{2}{c}{Values} \\
\cmidrule(l{3pt}r{3pt}){1-2} \cmidrule(l{3pt}r{3pt}){3-4}
T & N & s\_est\_1 & mise\_1 & hd\_1 & s\_0\_1 & s\_est\_2 & mise\_2 & hd\_2 & s\_0\_2 & s\_est\_3 & mise\_3 & hd\_3 & s\_0\_3\\
\midrule
33 & 25 & 0.426 & 0.2658961 & 0.3238095 & 0.580 & 0.694 & 0.3121780 & 8.951389 & 0.424 & 1.046 & 0.3426282 & 10.749311 & 0.274\\
33 & 50 & 0.752 & 0.1164306 & 0.0533333 & 0.250 & 1.522 & 0.1432683 & 4.298729 & 0.056 & 2.236 & 0.1705356 & 5.369388 & 0.020\\
33 & 100 & 0.998 & 0.0061372 & 0.0000000 & 0.002 & 1.984 & 0.0195518 & 0.176000 & 0.000 & 2.964 & 0.0364690 & 0.288000 & 0.000\\
33 & 200 & 1.000 & 0.0047217 & 0.0000000 & 0.000 & 2.000 & 0.0140943 & 0.000000 & 0.000 & 3.000 & 0.0275598 & 0.000000 & 0.000\\
65 & 25 & 0.364 & 0.2916927 & 0.9000000 & 0.640 & 0.630 & 0.3126324 & 16.958015 & 0.476 & 0.948 & 0.3495934 & 22.146628 & 0.318\\
\addlinespace
65 & 50 & 0.740 & 0.1182818 & 0.0704607 & 0.262 & 1.442 & 0.1529267 & 9.323210 & 0.078 & 2.112 & 0.1806322 & 11.904564 & 0.036\\
65 & 100 & 0.990 & 0.0060608 & 0.0000000 & 0.010 & 1.976 & 0.0114074 & 0.504000 & 0.000 & 2.956 & 0.0177055 & 0.704000 & 0.000\\
65 & 200 & 1.000 & 0.0013470 & 0.0000000 & 0.000 & 2.000 & 0.0039408 & 0.000000 & 0.000 & 3.000 & 0.0075983 & 0.000000 & 0.000\\
129 & 25 & 0.318 & 0.3168025 & 4.2565789 & 0.696 & 0.606 & 0.3226022 & 36.609375 & 0.488 & 0.858 & 0.3509816 & 48.355140 & 0.358\\
129 & 50 & 0.710 & 0.1307810 & 0.1977401 & 0.292 & 1.336 & 0.1734853 & 20.420091 & 0.124 & 2.014 & 0.1935275 & 26.887967 & 0.036\\
\addlinespace
129 & 100 & 0.968 & 0.0147474 & 0.0000000 & 0.032 & 1.966 & 0.0114233 & 1.462000 & 0.000 & 2.950 & 0.0133733 & 1.600000 & 0.000\\
129 & 200 & 1.000 & 0.0003877 & 0.0000000 & 0.000 & 2.000 & 0.0010812 & 0.000000 & 0.000 & 3.000 & 0.0021399 & 0.000000 & 0.000\\
\bottomrule
\end{tabular}