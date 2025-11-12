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
indices   = intersect(indices_0, indices_1)

# Initialize an empty results dataframe before the loop
all_comparison_results = data.frame()

if(file.exists('/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_read_id_sterile_1_trnd_0.rds')){
  inds_read = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_read_id_sterile_1_trnd_0.rds')
}else{
  inds_read = c()
}

inds2read = setdiff(indices,inds_read)

if(length(inds2read)>0){
  for (i in inds2read){
    message("Processing param_set_", i)
    results_0    = readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',i,'_sterile_1_trnd_0_tregs_0.rds'))
    results_1    = readRDS(paste0('./mass_sim_results_R/longitudinal_df_param_set_id_',i,'_sterile_1_trnd_0_tregs_1.rds'))
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
    min_reps  = min(full_data_comparison$rep_id)
    max_reps  = max(full_data_comparison$rep_id)
    t_max_ind = max(full_data_comparison$t)
    
    for (rep in min_reps:max_reps){  
      
      #### STERILE INJURY
      # tregs OFF
      full_data_comparison_scores_tregs_off = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on ==0)
      # tregs ON
      full_data_comparison_scores_tregs_on  = full_data_comparison %>% dplyr::filter(rep_id==rep & sterile==1 & tregs_on ==1)

      # --- Steady-state detection ---
      time_ss_tregs_off = steady_state_idx(full_data_comparison_scores_tregs_off$epithelial_score)
      time_ss_tregs_on  = steady_state_idx(full_data_comparison_scores_tregs_on$epithelial_score)
      time_ss_vec = c(time_ss_tregs_off, time_ss_tregs_on)
      
      if(!any(is.na(time_ss_vec))){
        # --- Comparisons --- without matching over steady state
        ## Treg OFF → ON 
        scores_tregs_off = full_data_comparison_scores_tregs_off$epithelial_score[time_ss_tregs_off:t_max_ind]
        scores_tregs_on  = full_data_comparison_scores_tregs_on$epithelial_score[time_ss_tregs_on:t_max_ind]
        
        cohens_d_on_off    = cohens_d(scores_tregs_off, scores_tregs_on)
        mean_on            = mean(scores_tregs_on)
        mean_off           = mean(scores_tregs_off)
        mean_diff_on_off   = mean_on - mean_off
        
        # --- Tabulate all comparisons ---
        comparison_results = data.frame(
          param_set_id = i,
          replicate_id = rep,
          comparison = c(
            "Treg_OFF_ON"
          ),
          injury_type = c("sterile"),
          ss_start_tregs_off = c(time_ss_tregs_off),
          ss_start_tregs_on  = c(time_ss_tregs_on),
          mean_treg_on       = c(mean_on),
          mean_treg_off      = c(mean_off),
          mean_diff_on_off   = c(mean_diff_on_off),
          cohens_d           = c(cohens_d_on_off)
        )
        
        # Append to global results
        all_comparison_results = bind_rows(all_comparison_results, comparison_results)
      }
    }
  }
  message("num of param_id successfully added: ", length(inds2read))
  
  #---- read previous file
  if(!file.exists('/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_read_id_sterile_1_trnd_0.rds')){
    saveRDS(inds2read, '/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_read_id_sterile_1_trnd_0.rds')
  }else{
    inds2read_old = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_read_id_sterile_1_trnd_0.rds')
    inds2read     = c(inds2read_old, inds2read)
    saveRDS(inds2read, '/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_read_id_sterile_1_trnd_0.rds')
  }
  
  #---- read previous file
  if(!file.exists('/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_sterile_1_trnd_0.rds')){
    saveRDS(all_comparison_results, '/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_sterile_1_trnd_0.rds')
  }else{
    all_comparison_results_old = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_sterile_1_trnd_0.rds')
    all_comparison_results     = rbind(all_comparison_results_old, all_comparison_results)
    saveRDS(all_comparison_results, '/Users/burcutepekule/Desktop/tregs/all_comparison_results_A13_sterile_1_trnd_0.rds')
  }
}else{
  message("No new pts added.")
}



