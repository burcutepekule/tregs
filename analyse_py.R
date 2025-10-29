library(tidyverse)

# Load all files
params    = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv')
results   = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/all_simulation_results.csv')

# Merge
full_data = results %>% left_join(params, by = 'param_set_id') 
