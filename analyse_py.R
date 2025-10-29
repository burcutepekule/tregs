library(tidyverse)

# Load all files
params    = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv')
results   = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/simulation_results_param_set_0.csv')

# Merge
full_data = results %>% left_join(params, by = 'param_set_id') 
