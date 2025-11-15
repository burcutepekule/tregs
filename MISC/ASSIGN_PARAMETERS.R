# ============================================================================
# ASSIGN PARAMETERS FROM CSV
# ============================================================================
# Thresholds
th_ROS_microbe = param_set_use$th_ROS_microbe
th_ROS_epith_recover = param_set_use$th_ROS_epith_recover
epith_recovery_chance = param_set_use$epith_recovery_chance
rat_com_pat_threshold = param_set_use$rat_com_pat_threshold

# Diffusion speeds
diffusion_speed_DAMPs = param_set_use$diffusion_speed_DAMPs
diffusion_speed_SAMPs = param_set_use$diffusion_speed_SAMPs
diffusion_speed_ROS = param_set_use$diffusion_speed_ROS

# Signal production
add_ROS = param_set_use$add_ROS
add_DAMPs = param_set_use$add_DAMPs
add_SAMPs = param_set_use$add_SAMPs

# Decay rates
ros_decay = param_set_use$ros_decay
DAMPs_decay = param_set_use$DAMPs_decay
SAMPs_decay = param_set_use$SAMPs_decay

# Activation thresholds
activation_threshold_DAMPs = param_set_use$activation_threshold_DAMPs
activation_threshold_SAMPs = param_set_use$activation_threshold_SAMPs

# Engulfment activities
activity_engulf_M0_baseline = param_set_use$activity_engulf_M0_baseline
activity_engulf_M1_baseline = param_set_use$activity_engulf_M1_baseline
activity_engulf_M2_baseline = param_set_use$activity_engulf_M2_baseline
activity_engulf_max = 0.99

# ROS production activities
activity_ROS_M0_baseline = 0.00
activity_ROS_M1_baseline = param_set_use$activity_ROS_M1_baseline
activity_ROS_M2_baseline = 0.00
activity_ROS_max = 0.99

# Leak rates
rate_leak_commensal_injury = param_set_use$rate_leak_commensal_injury
rate_leak_pathogen_injury = ifelse(sterile == 1, 0.0, param_set_use$rate_leak_pathogen_injury)
rate_leak_commensal_baseline = param_set_use$rate_leak_commensal_baseline

# Phagocyte parameters
active_age_limit = as.integer(param_set_use$active_age_limit)
cc_phagocyte = 5
digestion_time = 1

# Treg parameters
treg_vicinity_effect = 1
treg_discrimination_efficiency = param_set_use$treg_discrimination_efficiency
# allow_tregs_to_suppress_cognate = FALSE

# ============================================================================
# INITIALIZE SIMULATION
# ============================================================================
injury_site = get_middle_percent(seq(1, grid_size), injury_percentage)
n_pathogens_lp = round(rate_leak_pathogen_injury * length(injury_site))

precision = 10 * (exp(5 * treg_discrimination_efficiency))

# Calculate step sizes for activity increases
activity_engulf_M1_step = (activity_engulf_max - activity_engulf_M1_baseline) / cc_phagocyte
activity_engulf_M2_step = (activity_engulf_max - activity_engulf_M2_baseline) / cc_phagocyte
activity_ROS_M1_step = (activity_ROS_max - activity_ROS_M1_baseline) / cc_phagocyte