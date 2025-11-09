# rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
scenarios      = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results_presampled/scenarios.csv', show_col_types = FALSE)

# Load all files
params           = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results_presampled/loaded_parameters.csv', show_col_types = FALSE)


param_id     = 109 #6, 7
rep_ind      = 6

results      = read_csv(paste0('/Users/burcutepekule/Desktop/tregs/mass_sim_results_presampled/simulation_results_param_set_',param_id,'.csv'), show_col_types = FALSE)
results_plot = results %>% left_join(params, by = 'param_set_id') 

colnames(results_plot)[3] = 'tregs_on'

results_plot = results_plot %>% dplyr::filter(replicate_id==rep_ind)
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



