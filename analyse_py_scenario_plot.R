rm(list=ls())
library(dplyr)
library(tidyr)
source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
scenarios      = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/scenarios.csv', show_col_types = FALSE)

# Load all files
params           = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)
df_short_merged_with_params = readRDS('/Users/burcutepekule/Desktop/tregs/df_short_merged_with_params.rds')
table(df_short_merged_with_params$class_description)

df_check = df_short_merged_with_params %>% dplyr::filter(class_description=='Sterile better, pathogenic worse')

param_set_t = params %>% dplyr::filter(param_set_id %in% df_check$param_set_id) 
param_set_c = params %>% dplyr::filter(!(param_set_id %in% df_check$param_set_id))

for (ind in 1:dim(df_check)[1]){
  i            = df_check$param_set_id[ind]
  print(i)
  results      = read_csv(paste0('/Users/burcutepekule/Desktop/tregs/mass_sim_results/simulation_results_param_set_',i,'.csv'), show_col_types = FALSE)
  results_plot = results %>% left_join(params, by = 'param_set_id') 
  
  colnames(results_plot)[3] = 'tregs_on'
  
  # results_plot = results_plot %>% dplyr::filter(replicate_id==5)
  variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
  data_long = results_plot %>%
    dplyr::select(t, sterile, tregs_on, randomize_tregs, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
    geom_line(alpha = 1, linewidth = 1) +
    facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
  print(p)
  
  browser()
}


