rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
# Load the sampled parameters for the mass simulation. 
parameters_df= read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)
# > head(parameters_df)
# # A tibble: 6 × 25
# param_set_id th_ROS_microbe th_ROS_epith_recover epith_recovery_chance rat_com_pat_threshold diffusion_speed_DAMPs diffusion_speed_SAMPs diffusion_speed_ROS
# <dbl>          <dbl>                <dbl>                 <dbl>                 <dbl>                 <dbl>                 <dbl>               <dbl>
#   1            0          0.276               0.969                 0.0966                 0.522               0.0201                0.00637              0.0689
# 2            1          0.181               0.0695                0.409                  0.995               0.0822                0.103                0.0674
# 3            2          0.284               0.203                 0.0967                 0.927               0.00554               0.108                0.0858
# 4            3          0.730               0.445                 0.856                  0.794               0.0362                0.00896              0.0732
# 5            4          0.209               0.821                 0.413                  0.852               0.118                 0.0965               0.0165
# 6            5          0.280               0.509                 0.0383                 0.573               0.0942                0.0540               0.118 
# # ℹ 17 more variables: add_ROS <dbl>, add_DAMPs <dbl>, add_SAMPs <dbl>, ros_decay <dbl>, DAMPs_decay <dbl>, SAMPs_decay <dbl>, activation_threshold_DAMPs <dbl>,
# #   activation_threshold_SAMPs <dbl>, activity_engulf_M0_baseline <dbl>, activity_engulf_M1_baseline <dbl>, activity_engulf_M2_baseline <dbl>,
# #   activity_ROS_M1_baseline <dbl>, rate_leak_commensal_injury <dbl>, rate_leak_pathogen_injury <dbl>, rate_leak_commensal_baseline <dbl>, active_age_limit <dbl>,
# #   treg_discrimination_efficiency <dbl>

param_set_id = 1776 # parameter set index
rep_ind      = 5 # replication index, each parameter set has 10 replications to account for stochasticity

# Load the simulation outcome for a given param_set_id
results      = read_csv(paste0('/Users/burcutepekule/Desktop/tregs/mass_sim_results/simulation_results_param_set_',param_set_id,'.csv'), show_col_types = FALSE)

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
results      = results %>% dplyr::select(t, sterile, allow_tregs_to_do_their_job, randomize_tregs, parameter_names, endpoints)

# > print(results, n=10)
# # A tibble: 30,000 × 34
# t sterile allow_tregs_to_do_their_job randomize_tregs th_ROS_microbe th_ROS_epith_recover epith_recovery_chance rat_com_pat_threshold diffusion_speed_DAMPs
# <dbl>   <dbl>                       <dbl>           <dbl>          <dbl>                <dbl>                 <dbl>                 <dbl>                 <dbl>
#   1     0       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 2     1       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 3     2       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 4     3       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 5     4       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 6     5       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 7     6       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 8     7       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 9     8       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# 10     9       0                           0               0          0.615                0.147                0.0849                 0.518               0.00919
# # ℹ 29,990 more rows
# # ℹ 25 more variables: diffusion_speed_SAMPs <dbl>, diffusion_speed_ROS <dbl>, add_ROS <dbl>, add_DAMPs <dbl>, add_SAMPs <dbl>, ros_decay <dbl>, DAMPs_decay <dbl>,
# #   SAMPs_decay <dbl>, activation_threshold_DAMPs <dbl>, activation_threshold_SAMPs <dbl>, activity_engulf_M0_baseline <dbl>, activity_engulf_M1_baseline <dbl>,
# #   activity_engulf_M2_baseline <dbl>, activity_ROS_M1_baseline <dbl>, rate_leak_commensal_injury <dbl>, rate_leak_pathogen_injury <dbl>,
# #   rate_leak_commensal_baseline <dbl>, active_age_limit <dbl>, treg_discrimination_efficiency <dbl>, epithelial_healthy <dbl>, epithelial_inj_1 <dbl>,
# #   epithelial_inj_2 <dbl>, epithelial_inj_3 <dbl>, epithelial_inj_4 <dbl>, epithelial_inj_5 <dbl>
# # ℹ Use `print(n = ...)` to see more rows