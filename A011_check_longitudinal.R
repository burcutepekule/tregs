rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")

for(param_set_id in 11:100){
  df_treg_0_keep    = readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',param_set_id,'_tregs_0.rds'))
  df_treg_1_keep    = readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',param_set_id,'_tregs_1.rds'))
  
  df_treg   = rbind(df_treg_0_keep, df_treg_1_keep)
  
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
  
  print(sum(df_treg_0_keep$epithelial_healthy-df_treg_1_keep$epithelial_healthy))
  browser()
  
  # for(rep_id in 0:9){
    
    # df_treg_0 = df_treg_0_keep %>% filter(rep_id == rep_id)
    # df_treg_1 = df_treg_1_keep %>% filter(rep_id == rep_id)
    # df_treg   = rbind(df_treg_0, df_treg_1)
    # 
    # variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
    # data_long = df_treg %>%
    #   dplyr::select(t, sterile, tregs_on, randomize_tregs, all_of(variables)) %>%
    #   pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
    # 
    # p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
    #   geom_line(alpha = 1, linewidth = 1) +
    #   facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
    #   scale_color_manual(values = agent_colors) +
    #   theme_minimal() +
    #   labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
    # print(p)
    # 
    # print(sum(df_treg_0$epithelial_healthy-df_treg_1$epithelial_healthy))
    # browser()
  # }
}
