rm(list=ls())

library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(grid)
library(av)


dir_name_data = './mass_sim_results_R'
dir.create(dir_name_data, showWarnings = FALSE)

# Default values
# sterile_vec          = c(0,1) # 0 = infection, 1 = sterile injury
sterile_vec          = c(1) # 0 = infection, 1 = sterile injury
allow_tregs_vec      = c(0,1) # Allow tregs to do their job
# randomize_tregs_vec  = c(0,1) # 0 = follow DAMPs, 1 = random movement
randomize_tregs      = 0 # 0 = follow DAMPs, 1 = random movement
use_synchronized_rng = TRUE  # Use synchronized random numbers for fair comparisons

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================
params_df = read.csv("./original_lhs_parameters.csv", stringsAsFactors = FALSE)

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================
t_max      = 1000
plot_on    = 0
gif_on     = 0
plot_every = 1
grid_size  = 25
n_phagocytes = round(grid_size * grid_size * 0.35)
n_tregs = round(grid_size * grid_size * 0.35)
n_commensals_lp = 20

injury_percentage = 60
max_level_injury = 5

max_cell_value_ROS = 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1

lim_ROS = max_cell_value_ROS
lim_DAMP = max_cell_value_DAMPs
lim_SAMP = max_cell_value_SAMPs

act_radius_ROS = 1
act_radius_treg = 1
act_radius_DAMPs = 1
act_radius_SAMPs = 1

# Logistic function parameters (for epithelial injury calculation)
k_in = 0.044
x0_in = 50
shift_by = 10

param_set_id_use   = 5091 # 5091/2
rep_vec            = 2
# random_stream_file = paste0("./random_streams/random_numbers_seed_",param_set_id_use,".csv")
# stream_in          = scan(random_stream_file, quiet = TRUE, skip = 1)
# stream_in_long     = c()
# for(k in 1:150){
#   stream_in_long = c(stream_in_long, stream_in)
# }
# print(length(stream_in_long))
# saveRDS(stream_in_long, 'stream_in_long_5091.rds')
stream_in_long = readRDS('stream_in_long_5091.rds')

param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
param_set_read = param_set_use %>%
  pivot_longer(
    cols = -param_set_id,          # all columns except param_set_id
    names_to = "variable_name",
    values_to = "value"
  ) %>% dplyr::select(-param_set_id)

param_set_use$activity_engulf_M2_baseline=0.265 #even this is ok
# param_set_use$activity_engulf_M2_baseline=param_set_use$activity_engulf_M1_baseline

longitudinal_df_keep = c()
print(paste0('Processing param set ',param_set_id_use,' 😱'))

for(sterile in sterile_vec){
  for(allow_tregs in allow_tregs_vec){
    # for (reps_in in 0:9){
    
    dir_name = paste0('./frames_',param_set_id_use,'_sterile_',sterile,'_tregs_',allow_tregs)
    dir.create(dir_name, showWarnings = FALSE)
    
    for (reps_in in rep_vec){
      print(c(sterile, allow_tregs, reps_in))
      # rm(list = setdiff(ls(),c("param_set_id_use","param_set_read",
      #                          "reps_in","stream_in_long","sterile_vec",
      #                          "sterile","allow_tregs","randomize_tregs","allow_tregs_vec",
      #                          "use_synchronized_rng","params_df","param_set_use",
      #                          "randomize_tregs_vec","rep_vec","gif_on",
      #                          "colnames_insert","longitudinal_df_keep",
      #                          "dir_name_data","dir_name","t_max","plot_on",
      #                          "plot_every","grid_size","n_phagocytes","n_tregs",
      #                          "n_commensals_lp","injury_percentage","max_level_injury",
      #                          "max_cell_value_ROS","max_cell_value_DAMPs","max_cell_value_SAMPs",
      #                          "lim_ROS","lim_DAMP","lim_SAMP","act_radius_ROS","act_radius_treg",
      #                          "act_radius_DAMPs","act_radius_SAMPs","k_in","x0_in","shift_by")))
      
      # ============================================================================
      # LOAD FUNCTIONS
      # ============================================================================
      
      source("./MISC/FAST_FUNCTIONS.R")
      source("./MISC/PLOT_FUNCTIONS.R")
      source("./MISC/RUN_SIM_IN_A012.R")
      
      # ============================================================================
      # CREATE VIDEO
      # ============================================================================
      pattern = paste0("^frame_param_", param_set_id_use, "_rep_", reps_in,"_sterile_", sterile, "_tregs_", allow_tregs,
                       "_trnd_", randomize_tregs, "_\\d+\\.png$")
      
      png_files = list.files(dir_name, full.names = TRUE, pattern = pattern)
      png_files = png_files[order(as.numeric(gsub(".*_(\\d+)\\.png$", "\\1", png_files)))]
      video_out = paste0(dir_name, "/simulation_sterile", sterile, "_tregs_", allow_tregs,
                         "_trnd_", randomize_tregs, "_paramset_", param_set_id_use, ".mp4")
      
      if (length(png_files) > 0 & gif_on==1) {
        av_encode_video(
          input = png_files,
          output = video_out,
          framerate = 5,
          vfilter = "scale=1000:-2",
          codec = "libx264"
        )
        cat(sprintf("\nVideo created: %s\n", video_out))
      } else {
        cat("\nNo frames found for video creation.\n")
      }
    }
  }
}
colnames(longitudinal_df_keep)[c(7:37)] = colnames_insert

longitudinal_df_5000 = longitudinal_df_keep

df_treg   = longitudinal_df_5000 
variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
# variables = c('commensal','pathogen')
# variables = c(paste0("phagocyte_M1_L_", 0:5))
# variables = c(paste0("phagocyte_M2_L_", 0:5))

data_long = df_treg %>% dplyr::filter(sterile==1) %>%
  dplyr::select(t, sterile, tregs_on, randomize_tregs, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
  geom_line(alpha = 1, linewidth = 1) +
  facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
print(p)

