rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)

# Function to calculate Cohen's d
cohens_d = function(x, y) {
  nx = length(x)
  ny = length(y)
  mx = mean(x, na.rm = TRUE)
  my = mean(y, na.rm = TRUE)
  sx = sd(x, na.rm = TRUE)
  sy = sd(y, na.rm = TRUE)
  
  # Pooled standard deviation
  pooled_sd = sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  
  # Cohen's d
  d = (my - mx) / pooled_sd
  
  if (is.na(pooled_sd) || pooled_sd == 0) {
    return(0)  # or NA, depending on how you want to interpret it
  }
  return(d)
}

steady_state_idx = function(x, k = 20, tail_frac = 0.25,
                            tol_abs = 0.05*(150-25),     # 2.5 units
                            tol_sd  = 0.02*(150-25),# 1.25 units
                            tol_slope = 0.005*(150-25))  # 0.125/step
{
  n      = length(x)
  tail_n = ceiling(n*tail_frac)
  x_asym = mean(tail(x, tail_n))
  
  m    = rollapply(x, k, mean, align = "right", fill = NA)
  sd   = rollapply(x, k, sd,   align = "right", fill = NA)
  sl   = rollapply(x, k, function(v) mean(diff(v)), align = "right", fill = NA)
  cand = which(abs(m-x_asym) <= tol_abs & sd <= tol_sd & abs(sl) <= tol_slope)
  
  if (length(cand) == 0) return(NA_integer_)
  cand[1]
}

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
scenarios      = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/scenarios.csv', show_col_types = FALSE)

# Load all files
params                     = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)
results_merged             = c()
sterile_comparison_keep    = c()
pathogenic_comparison_keep = c()


path      = "/Users/burcutepekule/Desktop/tregs/mass_sim_results/"
files     = list.files(path, pattern = "^simulation_results_param_set_\\d+\\.csv$", full.names = TRUE)
indices   = str_extract(basename(files), "\\d+") |> as.numeric()
max_index = max(indices, na.rm = TRUE)
t_max_ind = 500 # max 500 time points

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

if(file.exists('/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id.rds')){
  start_id = 1+readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id.rds')
}else{
  start_id = 0
}

if(start_id<max_index){
  for (i in start_id:max_index){
    message("Processing param_set_", i)
    results         = read_csv(paste0('/Users/burcutepekule/Desktop/tregs/mass_sim_results/simulation_results_param_set_',i,'.csv'), show_col_types = FALSE)
    # results_merged  = rbind(results_merged, results)
    
    # Merge
    results = results %>% left_join(params, by = 'param_set_id') 
    results = results %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+ # higher the score, healthier the epithelium!
                                          5*epithelial_inj_1+
                                          4*epithelial_inj_2+
                                          3*epithelial_inj_3+
                                          2*epithelial_inj_4+
                                          1*epithelial_inj_5)
    
    full_data_comparison = results %>% dplyr::select(scenario_id, param_set_id, replicate_id, t, epithelial_score)
    
    for (rep in 0:99){ # PAY ATTENTION TO THIS!
      
      #### PATHOGENIC INJURY
      # tregs OFF
      full_data_comparison_scores_0 = full_data_comparison %>% dplyr::filter(replicate_id==rep & scenario_id==0)
      # tregs ON
      full_data_comparison_scores_1 = full_data_comparison %>% dplyr::filter(replicate_id==rep & scenario_id==1)
      # tregs ON, BUT ARE random
      full_data_comparison_scores_2 = full_data_comparison %>% dplyr::filter(replicate_id==rep & scenario_id==2)
      
      #### STERILE INJURY
      # tregs OFF
      full_data_comparison_scores_3 = full_data_comparison %>% dplyr::filter(replicate_id==rep & scenario_id==3)
      # tregs ON
      full_data_comparison_scores_4 = full_data_comparison %>% dplyr::filter(replicate_id==rep & scenario_id==4)
      # tregs ON, BUT ARE random
      full_data_comparison_scores_5 = full_data_comparison %>% dplyr::filter(replicate_id==rep & scenario_id==5)
      
      # --- Steady-state detection ---
      time_ss_0 = steady_state_idx(full_data_comparison_scores_0$epithelial_score)
      time_ss_1 = steady_state_idx(full_data_comparison_scores_1$epithelial_score)
      time_ss_2 = steady_state_idx(full_data_comparison_scores_2$epithelial_score)
      time_ss_3 = steady_state_idx(full_data_comparison_scores_3$epithelial_score)
      time_ss_4 = steady_state_idx(full_data_comparison_scores_4$epithelial_score)
      time_ss_5 = steady_state_idx(full_data_comparison_scores_5$epithelial_score)
      
      time_ss_vec = c(time_ss_0, time_ss_1, time_ss_2, time_ss_3, time_ss_4, time_ss_5)
      
      if(!any(is.na(time_ss_vec))){
        # --- Paired steady-state alignment points ---
        time_ss_01 = max(c(time_ss_0, time_ss_1)) # Treg OFF → ON
        time_ss_12 = max(c(time_ss_1, time_ss_2)) # Treg NOT RANDOM → RANDOM
        time_ss_34 = max(c(time_ss_3, time_ss_4)) # Treg OFF → ON
        time_ss_45 = max(c(time_ss_4, time_ss_5)) # Treg NOT RANDOM → RANDOM
        
        # --- Comparisons ---
        
        ## Treg OFF → ON (0 → 1)
        scores_01_0    = full_data_comparison_scores_0$epithelial_score[time_ss_01:t_max_ind]
        scores_01_1    = full_data_comparison_scores_1$epithelial_score[time_ss_01:t_max_ind]
        d_01           = cohens_d(scores_01_0, scores_01_1)
        mean_diff_01   = mean(scores_01_1) - mean(scores_01_0)
        time_vec       = time_ss_01:t_max_ind
        integrated_diff_01 = sum(diff(time_vec) * zoo::rollmean(scores_01_1 - scores_01_0, 2)) # integrate using trapezoidal rule
        
        ## Treg NOT RANDOM → RANDOM (1 → 2)
        scores_12_1    = full_data_comparison_scores_1$epithelial_score[time_ss_12:t_max_ind]
        scores_12_2    = full_data_comparison_scores_2$epithelial_score[time_ss_12:t_max_ind]
        d_12           = cohens_d(scores_12_1, scores_12_2)
        mean_diff_12   = mean(scores_12_2) - mean(scores_12_1)
        time_vec       = time_ss_12:t_max_ind
        integrated_diff_12 = sum(diff(time_vec) * zoo::rollmean(scores_12_2 - scores_12_1, 2)) # integrate using trapezoidal rule
        
        ## Treg OFF → ON (3 → 4)
        scores_34_3    = full_data_comparison_scores_3$epithelial_score[time_ss_34:t_max_ind]
        scores_34_4    = full_data_comparison_scores_4$epithelial_score[time_ss_34:t_max_ind]
        d_34           = cohens_d(scores_34_3, scores_34_4)
        mean_diff_34   = mean(scores_34_4) - mean(scores_34_3)
        time_vec       = time_ss_34:t_max_ind
        integrated_diff_34 = sum(diff(time_vec) * zoo::rollmean(scores_34_4 - scores_34_3, 2)) # integrate using trapezoidal rule
        
        ## Treg NOT RANDOM → RANDOM (4 → 5)
        scores_45_4    = full_data_comparison_scores_4$epithelial_score[time_ss_45:t_max_ind]
        scores_45_5    = full_data_comparison_scores_5$epithelial_score[time_ss_45:t_max_ind]
        d_45           = cohens_d(scores_45_4, scores_45_5)
        mean_diff_45   = mean(scores_45_5) - mean(scores_45_4)
        time_vec       = time_ss_45:t_max_ind
        integrated_diff_45 = sum(diff(time_vec) * zoo::rollmean(scores_45_5 - scores_45_4, 2)) # integrate using trapezoidal rule
        
        # --- Tabulate all comparisons ---
        comparison_results = data.frame(
          param_set_id = i,
          replicate_id = rep,
          comparison = c(
            "Treg_OFF_ON",
            "Treg_NRND_RND",
            "Treg_OFF_ON",
            "Treg_NRND_RND"
          ),
          injury_type = c("pathogenic", "pathogenic", "sterile", "sterile"),
          ss_start = c(time_ss_01, time_ss_12, time_ss_34, time_ss_45),
          cohens_d = c(d_01, d_12, d_34, d_45),
          mean_diff = c(mean_diff_01, mean_diff_12, mean_diff_34, mean_diff_45),
          integ_diff = c(integrated_diff_01, integrated_diff_12, integrated_diff_34, integrated_diff_45)
        )
        
        # Append to global results
        all_comparison_results = bind_rows(all_comparison_results, comparison_results)
      }
    }
  }
  message("Last param_id successfully added: ", max(all_comparison_results$param_set_id))
  
  saveRDS(max(all_comparison_results$param_set_id), '/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id.rds') 
  
  #---- read previous file
  if(!file.exists('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0.rds')){
    saveRDS(all_comparison_results, '/Users/burcutepekule/Desktop/tregs/all_comparison_results_0.rds')
  }else{
    all_comparison_results_old = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_0.rds')
    all_comparison_results = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results, '/Users/burcutepekule/Desktop/tregs/all_comparison_results_0.rds')
  }
}else{
  message("No new pts added.")
}



