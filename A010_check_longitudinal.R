rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

for(seed_in in 0:9){
  
  df_treg_0 = readRDS(paste0('longitudinal_df_',seed_in,'_tregs_0.rds'))
  df_treg_1 = readRDS(paste0('longitudinal_df_',seed_in,'_tregs_1.rds'))
  
  colnames(df_treg_0)[c(3,7:37)] = c('tregs_on',colnames_insert)
  colnames(df_treg_1)[c(3,7:37)] = c('tregs_on',colnames_insert)
  
  df_treg = rbind(df_treg_0, df_treg_1)
  
  variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
  data_long = df_treg %>%
    dplyr::select(t, sterile, tregs_on, randomize_tregs, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
    geom_line(alpha = 1, linewidth = 1) +
    facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
  print(p)
  
  print(sum(df_treg_0$epithelial_healthy-df_treg_1$epithelial_healthy))
  browser()
}


