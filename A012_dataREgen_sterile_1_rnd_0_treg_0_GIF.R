rm(list=ls())

library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(grid)
library(av)

dir_name = './frames'
dir.create(dir_name, showWarnings = FALSE)

dir_name_data = './mass_sim_results_R'
dir.create(dir_name_data, showWarnings = FALSE)

# Default values
sterile              = 1 # 0 = infection, 1 = sterile injury
allow_tregs_vec      = c(0,1)   # Allow tregs to do their job
randomize_tregs      = 0   # 0 = follow DAMPs, 1 = random movement
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
t_max      = 5000
plot_on    = 1
plot_every = 25
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
shift_by = 100

param_set_id_use = 5091
random_stream_file = paste0("./random_streams/random_numbers_seed_",param_set_id_use,".csv")
stream_in          = scan(random_stream_file, quiet = TRUE, skip = 1)
stream_in_long     = c()
for(k in 1:150){
  stream_in_long = c(stream_in_long, stream_in)
}
print(length(stream_in_long))


param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
longitudinal_df_keep = c()
print(paste0('Processing param set ',param_set_id_use,' 😱...'))

for(allow_tregs in allow_tregs_vec){
  # for (reps_in in 0:9){
  for (reps_in in 9){
    print(c(allow_tregs, reps_in))
    rm(list = setdiff(ls(),c("param_set_id_use","reps_in","stream_in",
                             "sterile","allow_tregs","randomize_tregs",
                             "use_synchronized_rng","params_df","param_set_use",
                             "colnames_insert","longitudinal_df_keep",
                             "dir_name_data","dir_name","t_max","plot_on",
                             "plot_every","grid_size","n_phagocytes","n_tregs",
                             "n_commensals_lp","injury_percentage","max_level_injury",
                             "max_cell_value_ROS","max_cell_value_DAMPs","max_cell_value_SAMPs",
                             "lim_ROS","lim_DAMP","lim_SAMP","act_radius_ROS","act_radius_treg",
                             "act_radius_DAMPs","act_radius_SAMPs","k_in","x0_in","shift_by")))
    
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

    if (length(png_files) > 0) {
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

colnames(longitudinal_df_keep)[c(7:37)] = colnames_insert

# saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
#                                      '_sterile_',sterile,
#                                      '_trnd_',randomize_tregs,
#                                      '_tregs_',allow_tregs,'.rds'))
# print(paste0('Data for param set id ',param_set_id_use,' saved 🥳.'))
longitudinal_df_500  = longitudinal_df_keep %>% dplyr::filter(t<=500)
longitudinal_df_1500 = longitudinal_df_keep %>% dplyr::filter(t<=1500)
longitudinal_df_5000 = longitudinal_df_keep

df_treg   = longitudinal_df_5000 %>% dplyr::filter(rep_id==9)
variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
variables = 'commensal'
variables = c(paste0("phagocyte_M1_L_", 0:5))
variables = c(paste0("phagocyte_M2_L_", 0:5))

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

######### CONFIRM
source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# ------------
results = longitudinal_df_500
results = results %>% left_join(params_df, by = 'param_set_id')
results = results %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+ # higher the score, healthier the epithelium!
                                      5*epithelial_inj_1+
                                      4*epithelial_inj_2+
                                      3*epithelial_inj_3+
                                      2*epithelial_inj_4+
                                      1*epithelial_inj_5)

full_data_comparison = results %>% dplyr::select(param_set_id, sterile, tregs_on, randomize_tregs, rep_id, t, epithelial_score)
min_reps  = min(full_data_comparison$rep_id)
max_reps  = max(full_data_comparison$rep_id)
t_max_ind = max(full_data_comparison$t)
all_comparison_results = c()
for (rep in min_reps:max_reps){

  #### STERILE INJURY
  # tregs OFF
  full_data_comparison_scores_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on ==0)
  # tregs ON
  full_data_comparison_scores_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on ==1)

  # --- Steady-state detection ---
  time_ss_0   = steady_state_idx(full_data_comparison_scores_0$epithelial_score)
  time_ss_1   = steady_state_idx(full_data_comparison_scores_1$epithelial_score)
  time_ss_vec = c(time_ss_0, time_ss_1)

  if(!any(is.na(time_ss_vec))){
    # --- Paired steady-state alignment points ---
    time_ss_01 = max(c(time_ss_0, time_ss_1)) # Treg OFF → ON

    # --- Comparisons ---
    ## Treg OFF → ON (3 → 4)
    scores_01_0    = full_data_comparison_scores_0$epithelial_score[time_ss_01:t_max_ind]
    scores_01_1    = full_data_comparison_scores_1$epithelial_score[time_ss_01:t_max_ind]
    d_01           = cohens_d(scores_01_0, scores_01_1)
    mean_diff_01   = mean(scores_01_1) - mean(scores_01_0)
    time_vec       = time_ss_01:t_max_ind
    integrated_diff_01 = sum(diff(time_vec) * zoo::rollmean(scores_01_1 - scores_01_0, 2)) # integrate using trapezoidal rule

    # --- Tabulate all comparisons ---
    comparison_results = data.frame(
      param_set_id = param_set_id_use,
      replicate_id = rep,
      comparison = c(
        "Treg_OFF_ON"
      ),
      injury_type = c("sterile"),
      ss_start    = c(time_ss_01),
      cohens_d    = c(d_01),
      mean_diff   = c(mean_diff_01),
      integ_diff  = c(integrated_diff_01)
    )

    # Append to global results
    all_comparison_results = bind_rows(all_comparison_results, comparison_results)
  }
}
all_comparison_results_500 = all_comparison_results

# ------------
results = longitudinal_df_1500
results = results %>% left_join(params_df, by = 'param_set_id')
results = results %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+ # higher the score, healthier the epithelium!
                                      5*epithelial_inj_1+
                                      4*epithelial_inj_2+
                                      3*epithelial_inj_3+
                                      2*epithelial_inj_4+
                                      1*epithelial_inj_5)

full_data_comparison = results %>% dplyr::select(param_set_id, sterile, tregs_on, randomize_tregs, rep_id, t, epithelial_score)
min_reps  = min(full_data_comparison$rep_id)
max_reps  = max(full_data_comparison$rep_id)
t_max_ind = max(full_data_comparison$t)
all_comparison_results = c()
for (rep in min_reps:max_reps){

  #### STERILE INJURY
  # tregs OFF
  full_data_comparison_scores_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on ==0)
  # tregs ON
  full_data_comparison_scores_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on ==1)

  # --- Steady-state detection ---
  time_ss_0   = steady_state_idx(full_data_comparison_scores_0$epithelial_score)
  time_ss_1   = steady_state_idx(full_data_comparison_scores_1$epithelial_score)
  time_ss_vec = c(time_ss_0, time_ss_1)

  if(!any(is.na(time_ss_vec))){
    # --- Paired steady-state alignment points ---
    time_ss_01 = max(c(time_ss_0, time_ss_1)) # Treg OFF → ON

    # --- Comparisons ---
    ## Treg OFF → ON (3 → 4)
    scores_01_0    = full_data_comparison_scores_0$epithelial_score[time_ss_01:t_max_ind]
    scores_01_1    = full_data_comparison_scores_1$epithelial_score[time_ss_01:t_max_ind]
    d_01           = cohens_d(scores_01_0, scores_01_1)
    mean_diff_01   = mean(scores_01_1) - mean(scores_01_0)
    time_vec       = time_ss_01:t_max_ind
    integrated_diff_01 = sum(diff(time_vec) * zoo::rollmean(scores_01_1 - scores_01_0, 2)) # integrate using trapezoidal rule

    # --- Tabulate all comparisons ---
    comparison_results = data.frame(
      param_set_id = param_set_id_use,
      replicate_id = rep,
      comparison = c(
        "Treg_OFF_ON"
      ),
      injury_type = c("sterile"),
      ss_start    = c(time_ss_01),
      cohens_d    = c(d_01),
      mean_diff   = c(mean_diff_01),
      integ_diff  = c(integrated_diff_01)
    )

    # Append to global results
    all_comparison_results = bind_rows(all_comparison_results, comparison_results)
  }
}
all_comparison_results_1500 = all_comparison_results


# ------------
results = longitudinal_df_5000
results = results %>% left_join(params_df, by = 'param_set_id')
results = results %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+ # higher the score, healthier the epithelium!
                                      5*epithelial_inj_1+
                                      4*epithelial_inj_2+
                                      3*epithelial_inj_3+
                                      2*epithelial_inj_4+
                                      1*epithelial_inj_5)

full_data_comparison = results %>% dplyr::select(param_set_id, sterile, tregs_on, randomize_tregs, rep_id, t, epithelial_score)
min_reps  = min(full_data_comparison$rep_id)
max_reps  = max(full_data_comparison$rep_id)
t_max_ind = max(full_data_comparison$t)
all_comparison_results = c()
for (rep in min_reps:max_reps){

  #### STERILE INJURY
  # tregs OFF
  full_data_comparison_scores_0 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on ==0)
  # tregs ON
  full_data_comparison_scores_1 = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on ==1)

  # --- Steady-state detection ---
  time_ss_0   = steady_state_idx(full_data_comparison_scores_0$epithelial_score)
  time_ss_1   = steady_state_idx(full_data_comparison_scores_1$epithelial_score)
  time_ss_vec = c(time_ss_0, time_ss_1)

  if(!any(is.na(time_ss_vec))){
    # --- Paired steady-state alignment points ---
    time_ss_01 = max(c(time_ss_0, time_ss_1)) # Treg OFF → ON

    # --- Comparisons ---
    ## Treg OFF → ON (3 → 4)
    scores_01_0    = full_data_comparison_scores_0$epithelial_score[time_ss_01:t_max_ind]
    scores_01_1    = full_data_comparison_scores_1$epithelial_score[time_ss_01:t_max_ind]
    d_01           = cohens_d(scores_01_0, scores_01_1)
    mean_diff_01   = mean(scores_01_1) - mean(scores_01_0)
    time_vec       = time_ss_01:t_max_ind
    integrated_diff_01 = sum(diff(time_vec) * zoo::rollmean(scores_01_1 - scores_01_0, 2)) # integrate using trapezoidal rule

    # --- Tabulate all comparisons ---
    comparison_results = data.frame(
      param_set_id = param_set_id_use,
      replicate_id = rep,
      comparison = c(
        "Treg_OFF_ON"
      ),
      injury_type = c("sterile"),
      ss_start    = c(time_ss_01),
      cohens_d    = c(d_01),
      mean_diff   = c(mean_diff_01),
      integ_diff  = c(integrated_diff_01)
    )

    # Append to global results
    all_comparison_results = bind_rows(all_comparison_results, comparison_results)
  }
}
all_comparison_results_5000 = all_comparison_results