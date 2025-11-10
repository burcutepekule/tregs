rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)


source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# Load files
params                     = read_csv('/Users/burcutepekule/Desktop/tregs/original_lhs_parameters.csv', show_col_types = FALSE)
results_merged             = c()
sterile_comparison_keep    = c()
pathogenic_comparison_keep = c()


path      = "/Users/burcutepekule/Desktop/tregs/mass_sim_results_R/"
files_0   = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_1_trnd_0_tregs_0.rds$", full.names = TRUE)
files_1   = list.files(path, pattern = "^longitudinal_df_param_set_id_\\d+\\_sterile_1_trnd_0_tregs_1.rds$", full.names = TRUE)
indices_0 = str_extract(basename(files_0), "\\d+") |> as.numeric()
indices_1 = str_extract(basename(files_1), "\\d+") |> as.numeric()
max_index = min(max(indices_0), max(indices_1))
t_max_ind = 500 # max 500 time points

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

if(file.exists('/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id_sterile_1_trnd_0.rds')){
  start_id = 1+readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id_sterile_1_trnd_0.rds')
}else{
  start_id = 0
}

if(start_id<max_index){
  for (i in start_id:max_index){
    message("Processing param_set_", i)
    results_0    = readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',i,'_tregs_0.rds'))
    results_1    = readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',i,'_tregs_1.rds'))
    results      = rbind(results_0, results_1)

    # Merge
    results = results %>% left_join(params, by = 'param_set_id') 
    results = results %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+ # higher the score, healthier the epithelium!
                                          5*epithelial_inj_1+
                                          4*epithelial_inj_2+
                                          3*epithelial_inj_3+
                                          2*epithelial_inj_4+
                                          1*epithelial_inj_5)
    
    full_data_comparison = results %>% dplyr::select(param_set_id, sterile, tregs_on, randomize_tregs, rep_id, t, epithelial_score)
    min_reps = min(full_data_comparison$rep_id)
    max_reps = max(full_data_comparison$rep_id)
    
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
          param_set_id = i,
          replicate_id = rep,
          comparison = c(
            "Treg_OFF_ON",
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
  }
  message("Last param_id successfully added: ", max(all_comparison_results$param_set_id))
  
  saveRDS(max(all_comparison_results$param_set_id), '/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id_sterile_1_trnd_0.rds') 
  
  #---- read previous file
  if(!file.exists('/Users/burcutepekule/Desktop/tregs/all_comparison_results_sterile_1_trnd_0.rds')){
    saveRDS(all_comparison_results, '/Users/burcutepekule/Desktop/tregs/all_comparison_results_sterile_1_trnd_0.rds')
  }else{
    all_comparison_results_old = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_sterile_1_trnd_0.rds')
    all_comparison_results = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results, '/Users/burcutepekule/Desktop/tregs/all_comparison_results_sterile_1_trnd_0.rds')
  }
}else{
  message("No new pts added.")
}



