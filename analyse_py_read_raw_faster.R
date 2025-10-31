rm(list = ls())
library(dplyr)
library(readr)
library(stringr)
library(zoo)
library(data.table)   # = for faster I/O and filtering

# --- Functions ----------------------------------------------------------

cohens_d = function(x, y) {
  nx = length(x)
  ny = length(y)
  mx = mean(x, na.rm = TRUE)
  my = mean(y, na.rm = TRUE)
  sx = sd(x, na.rm = TRUE)
  sy = sd(y, na.rm = TRUE)
  pooled_sd = sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  if (is.na(pooled_sd) || pooled_sd == 0) return(0)
  (my - mx) / pooled_sd
}

steady_state_idx = function(x, k = 20, tail_frac = 0.1,
                            tol_abs = 0.02 * (150 - 25),
                            tol_sd  = 0.5  * 0.02 * (150 - 25),
                            tol_slope = 0.005 * (150 - 25)) {
  n = length(x)
  tail_n = max(ceiling(n * tail_frac), k)
  x_asym = median(tail(x, tail_n), na.rm = TRUE)
  m  = rollmean(x, k, align = "right", fill = NA)
  sd = rollapply(x, k, sd, align = "right", fill = NA)
  sl = rollapply(x, k, function(v) mean(diff(v)), align = "right", fill = NA)
  cand = which(abs(m - x_asym) <= tol_abs & sd <= tol_sd & abs(sl) <= tol_slope)
  if (length(cand) == 0) return(NA_integer_)
  cand[1]
}

# --- Setup --------------------------------------------------------------

params = fread('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv')
path   = "/Users/burcutepekule/Desktop/tregs/mass_sim_results/"
files  = list.files(path, pattern = "^simulation_results_param_set_\\d+\\.csv$", full.names = TRUE)
indices = as.numeric(str_extract(basename(files), "\\d+"))
max_index = max(indices, na.rm = TRUE)
t_max_ind = 500
t_cut = 0

all_comparison_results = list()  # preallocate list instead of growing a df

if(file.exists('/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id.rds')){
  start_id = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id.rds')
}else{
  start_id = 0
}

# --- Main Loop ----------------------------------------------------------
for (i in start_id:max_index){
  message("Processing param_set_", i)
  
  # fread is ~10x faster than read_csv
  results = fread(file.path(path, paste0("simulation_results_param_set_", i, ".csv")))
  results = merge(results, params, by = "param_set_id", all.x = TRUE)
  
  # Precompute epithelial_score directly in data.table
  results[, epithelial_score := 6 * epithelial_healthy +
            5 * epithelial_inj_1 +
            4 * epithelial_inj_2 +
            3 * epithelial_inj_3 +
            2 * epithelial_inj_4 +
            1 * epithelial_inj_5]
  
  # Select only necessary columns
  full_data_comparison = results[, .(scenario_id, param_set_id, replicate_id, t, epithelial_score)]
  
  for (rep in 0:9) {
    # Pre-slice all 6 scenario subsets once using data.table filtering (fast)
    subset_dt = full_data_comparison[replicate_id == rep]
    sc = split(subset_dt, subset_dt$scenario_id)
    
    # Skip if any scenario missing
    if (length(sc) < 6) next
    
    ss_times = sapply(0:5, function(s) steady_state_idx(sc[[as.character(s)]]$epithelial_score))
    if (any(is.na(ss_times))) next
    
    # Define paired steady-state points
    time_ss_01 = max(ss_times[c(1, 2)])
    time_ss_12 = max(ss_times[c(2, 3)])
    time_ss_34 = max(ss_times[c(4, 5)])
    time_ss_45 = max(ss_times[c(5, 6)])
    
    # Helper for pair comparisons
    get_metrics = function(a, b, t_ss) {
      x = a$epithelial_score[t_ss:t_max_ind]
      y = b$epithelial_score[t_ss:t_max_ind]
      list(
        d = cohens_d(x, y),
        med = median(y) - median(x),
        mean = mean(y) - mean(x)
      )
    }
    
    res_01 = get_metrics(sc[["0"]], sc[["1"]], time_ss_01)
    res_12 = get_metrics(sc[["1"]], sc[["2"]], time_ss_12)
    res_34 = get_metrics(sc[["3"]], sc[["4"]], time_ss_34)
    res_45 = get_metrics(sc[["4"]], sc[["5"]], time_ss_45)
    
    # Build small result table for this replicate
    comparison_results = data.frame(
      param_set_id = i,
      replicate_id = rep,
      comparison = c("Treg_OFF_ON", "Treg_NRND_RND", "Treg_OFF_ON", "Treg_NRND_RND"),
      injury_type = c("pathogenic", "pathogenic", "sterile", "sterile"),
      ss_start = c(time_ss_01, time_ss_12, time_ss_34, time_ss_45),
      cohens_d = c(res_01$d, res_12$d, res_34$d, res_45$d),
      median_diff = c(res_01$med, res_12$med, res_34$med, res_45$med),
      mean_diff = c(res_01$mean, res_12$mean, res_34$mean, res_45$mean)
    )
    
    all_comparison_results[[length(all_comparison_results) + 1]] = comparison_results
  }
}

# --- Combine and save ---------------------------------------------------

all_comparison_results = data.table::rbindlist(all_comparison_results)
saveRDS(all_comparison_results, paste0('/Users/burcutepekule/Desktop/tregs/all_comparison_results_fast_',start_id,'.rds'))
saveRDS(max(all_comparison_results$param_set_id), '/Users/burcutepekule/Desktop/tregs/all_comparison_results_max_id.rds')

