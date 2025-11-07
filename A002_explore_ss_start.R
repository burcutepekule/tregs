rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")

df_raw_keep = c()
for(ss_start_threshold in c(50,100,150,200,250,300,350,400,450)){
  df_raw       = readRDS(paste0('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_raw_',ss_start_threshold,'.rds'))
  df_raw$ss_th = ss_start_threshold
  df_raw_keep = rbind(df_raw_keep, df_raw)
}

plot(c(50,100,150,200,250,300,350,400,450), table(df_raw_keep$ss_th))

unique_to_450 = df_raw_keep %>%
  group_by(param_set_id, replicate_id) %>%
  summarise(
    ss_th_values = list(unique(ss_th)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  filter(length(ss_th_values) == 1 && all(ss_th_values == 450)) %>%
  ungroup() %>%
  select(param_set_id, replicate_id)

