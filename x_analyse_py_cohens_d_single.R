rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv

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
  return(d)
}

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
scenarios      = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/scenarios.csv', show_col_types = FALSE)

# Load all files
params                     = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)
results_merged             = c()
sterile_comparison_keep    = c()
pathogenic_comparison_keep = c()

t_cut = 0
i=22
print(i)
results         = read_csv(paste0('/Users/burcutepekule/Desktop/tregs/mass_sim_results/simulation_results_param_set_',i,'.csv'), show_col_types = FALSE)
results = results %>% left_join(params, by = 'param_set_id') 

# ----- PLOT TO CHECK FIRST ---------------------------------------------------------------
results_plot = results
colnames(results_plot)[3] = 'tregs_on'
# results_plot = results_plot %>% dplyr::filter(replicate_id==5)
variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
data_long = results_plot %>%
  dplyr::select(t, sterile, tregs_on, randomize_tregs, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
  geom_line(alpha = 1, linewidth = 1) +
  facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
print(p)
# ----- PLOT TO CHECK FIRST ---------------------------------------------------------------

# Merge
results = results %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+
                                      5*epithelial_inj_1+
                                      4*epithelial_inj_2+
                                      3*epithelial_inj_3+
                                      2*epithelial_inj_4+
                                      1*epithelial_inj_5)

# Keep ALL data points including all replicates for t>=249
full_data_comparison = results %>%
  filter(t >= t_cut) %>%
  dplyr::select(scenario_id, param_set_id, replicate_id, t, epithelial_score)

pathogenic_comparison = full_data_comparison %>%
  filter(scenario_id %in% c(0, 1)) %>%
  group_by(param_set_id, replicate_id) %>%
  nest() %>%
  mutate(
    test_result = map(data, function(df) {
      scores_0 = df %>% filter(scenario_id == 0) %>% pull(epithelial_score)
      scores_1 = df %>% filter(scenario_id == 1) %>% pull(epithelial_score)
      
      # Only run tests if we have data for both scenarios
      if(length(scores_0) > 0 && length(scores_1) > 0) {
        d = cohens_d(scores_0, scores_1)
        # Calculate effect size (median difference)
        median_diff = median(scores_1, na.rm = TRUE) - median(scores_0, na.rm = TRUE)
        mean_diff = mean(scores_1, na.rm = TRUE) - mean(scores_0, na.rm = TRUE)
        
        return(tibble(
          cohens_d = d,
          median_diff = median_diff,
          mean_diff = mean_diff,
          n_points_0 = length(scores_0),
          n_points_1 = length(scores_1)
        ))
      } else {
        return(tibble(
          cohens_d = NA_real_,
          median_diff = NA_real_,
          mean_diff = NA_real_,
          n_points_0 = length(scores_0),
          n_points_1 = length(scores_1)
        ))
      }
    })
  ) %>%
  unnest(test_result) %>%
  dplyr::select(-data) %>%
  mutate(condition = "Sterile",
         effect_size = case_when(
           abs(cohens_d) < 0.2 ~ "Negligible",
           abs(cohens_d) < 0.5 ~ "Small",
           abs(cohens_d) < 0.8 ~ "Medium",
           TRUE ~ "Large"
         )
  )

# For pathogenic condition (3 vs 4)
sterile_comparison = full_data_comparison %>%
  filter(scenario_id %in% c(3, 4)) %>%
  group_by(param_set_id, replicate_id) %>%
  nest() %>%
  mutate(
    test_result = map(data, function(df) {
      scores_3 = df %>% filter(scenario_id == 3) %>% pull(epithelial_score)
      scores_4 = df %>% filter(scenario_id == 4) %>% pull(epithelial_score)
      
      if(length(scores_3) > 0 && length(scores_4) > 0) {
        d = cohens_d(scores_3, scores_4)
        median_diff = median(scores_4, na.rm = TRUE) - median(scores_3, na.rm = TRUE)
        mean_diff = mean(scores_4, na.rm = TRUE) - mean(scores_3, na.rm = TRUE)
        
        return(tibble(
          cohens_d = d,
          median_diff = median_diff,
          mean_diff = mean_diff,
          n_points_3 = length(scores_3),
          n_points_4 = length(scores_4)
        ))
      } else {
        return(tibble(
          cohens_d = NA_real_,
          median_diff = NA_real_,
          mean_diff = NA_real_,
          n_points_3 = length(scores_3),
          n_points_4 = length(scores_4)
        ))
      }
    })
  ) %>%
  unnest(test_result) %>%
  dplyr::select(-data) %>%
  mutate(condition = "Pathogenic",
         effect_size = case_when(
           abs(cohens_d) < 0.2 ~ "Negligible",
           abs(cohens_d) < 0.5 ~ "Small",
           abs(cohens_d) < 0.8 ~ "Medium",
           TRUE ~ "Large"
         )
  )

sterile_comparison_keep    = rbind(sterile_comparison_keep, sterile_comparison)
pathogenic_comparison_keep = rbind(pathogenic_comparison_keep, pathogenic_comparison)

table(sterile_comparison_keep$effect_size)
table(pathogenic_comparison_keep$effect_size)

sterile_comparison_keep    = sterile_comparison_keep %>% dplyr::mutate(sterile_treg_better = ifelse(((effect_size=='Large' | effect_size=='Medium') & mean_diff>0 & median_diff>0),1,0))
pathogenic_comparison_keep = pathogenic_comparison_keep %>% dplyr::mutate(pathogenic_treg_better = ifelse(((effect_size=='Large' | effect_size=='Medium') & mean_diff>0 & median_diff>0),1,0))

sterile_comparison_keep    = sterile_comparison_keep %>% dplyr::mutate(sterile_treg_worse = ifelse(((effect_size=='Large' | effect_size=='Medium') & mean_diff<0 & median_diff<0),1,0))
pathogenic_comparison_keep = pathogenic_comparison_keep %>% dplyr::mutate(pathogenic_treg_worse = ifelse(((effect_size=='Large' | effect_size=='Medium') & mean_diff<0 & median_diff<0),1,0))

sterile_comparison_keep_short    = sterile_comparison_keep %>% dplyr::select(param_set_id, replicate_id, sterile_treg_better, sterile_treg_worse)
pathogenic_comparison_keep_short = pathogenic_comparison_keep %>% dplyr::select(param_set_id, replicate_id, pathogenic_treg_better, pathogenic_treg_worse)
df_short = merge(sterile_comparison_keep_short, pathogenic_comparison_keep_short, by=c('param_set_id','replicate_id'))

df_short = df_short %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::mutate(
    sterile_treg_better_avg = mean(sterile_treg_better, na.rm = TRUE),
    pathogenic_treg_better_avg = mean(pathogenic_treg_better, na.rm = TRUE),
    sterile_treg_worse_avg = mean(sterile_treg_worse, na.rm = TRUE),
    pathogenic_treg_worse_avg = mean(pathogenic_treg_worse, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

df_short = df_short %>% dplyr::select(param_set_id, sterile_treg_better_avg, pathogenic_treg_better_avg, sterile_treg_worse_avg, pathogenic_treg_worse_avg)
df_short = unique(df_short)

df_short_no_impact = df_short %>% dplyr::filter(sterile_treg_better_avg<=0.5 & pathogenic_treg_better_avg<=0.5 & sterile_treg_worse_avg<=0.5 & pathogenic_treg_worse_avg<=0.5)
df_short_impact    = df_short %>% dplyr::filter(!(sterile_treg_better_avg<=0.5 & pathogenic_treg_better_avg<=0.5 & sterile_treg_worse_avg<=0.5 & pathogenic_treg_worse_avg<=0.5))


df_short_impact = df_short_impact %>%
  mutate(
    class_code = paste0(
      as.numeric(sterile_treg_better_avg > 0.5),
      as.numeric(pathogenic_treg_better_avg > 0.5),
      as.numeric(sterile_treg_worse_avg > 0.5),
      as.numeric(pathogenic_treg_worse_avg > 0.5)
    ),
    class_description = case_when(
      class_code == "0000" ~ "No effect",
      class_code == "1000" ~ "Sterile better only",
      class_code == "0100" ~ "Pathogenic better only",
      class_code == "0010" ~ "Sterile worse only",
      class_code == "0001" ~ "Pathogenic worse only",
      class_code == "1100" ~ "Both better",
      class_code == "0011" ~ "Both worse",
      class_code == "1010" ~ "Sterile conflicting",
      class_code == "0101" ~ "Pathogenic conflicting",
      class_code == "1001" ~ "Sterile better, pathogenic worse",
      class_code == "0110" ~ "Pathogenic better, sterile worse",
      class_code == "1110" ~ "Pathogenic better, sterile conflicting",
      class_code == "1101" ~ "Sterile better, pathogenic conflicting",
      class_code == "1011" ~ "Sterile conflicting, pathogenic worse",
      class_code == "0111" ~ "Pathogenic conflicting, sterile worse",
      class_code == "1111" ~ "All conflicting",
      TRUE ~ "Other"
    )
  )

# align
df_short_no_impact$class_code = '0000'
df_short_no_impact$class_description ='No effect'

df_short_merged = rbind(df_short_no_impact, df_short_impact)
table(df_short_merged$class_description)

df_short_merged_with_params = merge(df_short_merged, params, by='param_set_id')



