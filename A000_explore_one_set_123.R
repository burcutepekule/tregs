rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
# Load the sampled parameters for the mass simulation. 
parameters_df= read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results_presampled_123/loaded_parameters.csv', show_col_types = FALSE)

param_set_id = 123 # parameter set index
rep_ind      = 61 # replication index, each parameter set has 10 replications to account for stochasticity
# compare rep_ind 43 va 61 -> very similar start_ss, very different outcomes


# Load the simulation outcome for a given param_set_id
results      = read_csv(paste0('/Users/burcutepekule/Desktop/tregs/mass_sim_results_presampled_123/simulation_results_param_set_',param_set_id,'.csv'), show_col_types = FALSE)
# Add the parameters to the dataframe by merhing over param_set_id
results      = results %>% left_join(parameters_df, by = 'param_set_id') 

# Pick only the outcome we are interested in, which is the epithelial cell types over time
# with 6 different health (or injury) levels (_healthy, _inj_1, _inj_2, _inj_3, _inj_4, _inj_5) 
endpoints    = c("epithelial_healthy", paste0("epithelial_inj_", 1:5)) # we are interested in these time series

# All parameter names
parameter_names = c(
  "th_ROS_microbe",
  "th_ROS_epith_recover",
  "epith_recovery_chance",
  "rat_com_pat_threshold",
  "diffusion_speed_DAMPs",
  "diffusion_speed_SAMPs",
  "diffusion_speed_ROS",
  "add_ROS",
  "add_DAMPs",
  "add_SAMPs",
  "ros_decay",
  "DAMPs_decay",
  "SAMPs_decay",
  "activation_threshold_DAMPs",
  "activation_threshold_SAMPs",
  "activity_engulf_M0_baseline",
  "activity_engulf_M1_baseline",
  "activity_engulf_M2_baseline",
  "activity_ROS_M1_baseline",
  "rate_leak_commensal_injury",
  "rate_leak_pathogen_injury",
  "rate_leak_commensal_baseline",
  "active_age_limit",
  "treg_discrimination_efficiency"
)

results_plot = results
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
