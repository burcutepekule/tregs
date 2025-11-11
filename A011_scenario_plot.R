rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")

# Load all files

param_id     = 5091
rep_ind      = 8 # 0, 1, 3, 4, 5, 9 Tregs better, 2, Tregs worse, 6, 7, 8 same

results_0    = readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',param_id,'_sterile_1_trnd_0_tregs_0.rds'))
results_1    = readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',param_id,'_sterile_1_trnd_0_tregs_1.rds'))
results      = rbind(results_0, results_1)

results = results %>% dplyr::filter(rep_id==rep_ind)
variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
data_long = results %>%
  dplyr::select(t, sterile, tregs_on, randomize_tregs, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
  geom_line(alpha = 1, linewidth = 1) +
  facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
print(p)

df_raw    = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_sterile_1_trnd_0.rds')
df_raw_pick = df_raw %>% dplyr::filter(param_set_id==param_id)

df_params = read_csv('/Users/burcutepekule/Desktop/tregs/original_lhs_parameters.csv', show_col_types = FALSE)
df_params_pick = df_params %>%
  filter(param_set_id == param_id) %>%
  pivot_longer(
    cols = -param_set_id,          # all columns except param_set_id
    names_to = "parameter",
    values_to = "value"
  )



