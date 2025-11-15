rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
inj_type             = 'sterile'
ss_start_threshold   = 4500
t_max                = 5000
tol_in               = 25*0.25

df_raw    = readRDS('/Users/burcutepekule/Desktop/tregs/all_comparison_results_A15_sterile_1_trnd_0.rds')
length(unique(df_raw$param_set_id))
df_params = read_csv('/Users/burcutepekule/Desktop/tregs/original_lhs_parameters.csv', show_col_types = FALSE)

df_plot = df_raw[c('param_set_id','replicate_id',
                    'ss_start_tregs_off','ss_start_tregs_on',
                    'mean_treg_on','mean_treg_off')]

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
param_id_all_below = df_plot %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start_tregs_off <= ss_start_threshold & ss_start_tregs_on <= ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_plot = df_plot %>% dplyr::filter(param_set_id %in% param_id_all_below)

df_summary = df_plot %>%
  group_by(param_set_id) %>%
  summarise(
    mean_treg_on_mean  = mean(mean_treg_on,  na.rm = TRUE),
    mean_treg_off_mean = mean(mean_treg_off, na.rm = TRUE),
    mean_treg_on_off_diff_mean = mean_treg_on_mean-mean_treg_off_mean
  )

tol_in = 25*0.25
df_summary = df_summary %>% dplyr::mutate(rounded_on_off_diff_mean = ifelse(
  mean_treg_on_off_diff_mean < tol_in & mean_treg_on_off_diff_mean > -1*tol_in, 
  0, mean_treg_on_off_diff_mean
))

param_ids_check = df_summary$param_set_id

df_raw_old    = readRDS('/Users/burcutepekule/Desktop/tregs/results_A13/all_comparison_results_A13_sterile_1_trnd_0.rds')

df_plot = df_raw_old[c('param_set_id','replicate_id',
                   'ss_start_tregs_off','ss_start_tregs_on',
                   'mean_treg_on','mean_treg_off')]

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
param_id_all_below = df_plot %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start_tregs_off <= ss_start_threshold & ss_start_tregs_on <= ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
df_plot = df_plot %>% dplyr::filter(param_set_id %in% param_id_all_below)

df_summary = df_plot %>%
  group_by(param_set_id) %>%
  summarise(
    mean_treg_on_mean  = mean(mean_treg_on,  na.rm = TRUE),
    mean_treg_off_mean = mean(mean_treg_off, na.rm = TRUE),
    mean_treg_on_off_diff_mean = mean_treg_on_mean-mean_treg_off_mean
  )

tol_in = 25*0.25
df_summary = df_summary %>% dplyr::mutate(rounded_on_off_diff_mean = ifelse(
  mean_treg_on_off_diff_mean < tol_in & mean_treg_on_off_diff_mean > -1*tol_in, 
  0, mean_treg_on_off_diff_mean
))

df_summary = df_summary %>% dplyr::filter(param_set_id %in% param_ids_check)


















hist(df_raw$mean_treg_off)
hist(df_raw$mean_treg_on)

df_raw_keep = df_raw
df_raw = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)
hist(df_raw_keep$mean_diff)



df_raw     = df_raw %>% dplyr::mutate(abs_cohens_d = abs(cohens_d))
df_raw     = df_raw %>% dplyr::mutate(effect_size = case_when(
  abs_cohens_d < 0.2 ~ "Negligible",
  abs_cohens_d < 0.5 & abs_cohens_d>= 0.2  ~ "Small",
  abs_cohens_d < 0.8 & abs_cohens_d>= 0.5 ~ "Medium",
  TRUE ~ "Large"
))

df_raw_plot_nz = df_raw %>% dplyr::filter(effect_size %in% c("Medium","Large") & (mean_diff_on_off>tol_in | mean_diff_on_off<(-1*tol_in)))
df_raw_plot_z  = df_raw %>% dplyr::filter((mean_diff_on_off<=tol_in & mean_diff_on_off>=(-1*tol_in))) # this is not true- what about effect_size small but abs(mean_diff)>tol_in?
df_raw_plot    = rbind(df_raw_plot_z, df_raw_plot_nz)
hist(df_raw_plot$mean_diff,30)

x   = df_raw_plot$mean_diff
round(100*sum(x>tol_in)/length(x),2)
round(100*sum(x<=tol_in & x>=(-1*tol_in))/length(x),2)
round(100*sum(x<(-1*tol_in))/length(x),2)
hist(x,30)

#-----------
# var_df = df_raw_plot %>%
var_df = df_raw %>%
  group_by(param_set_id) %>%
  summarise(
    sd_on = sd(mean_treg_on),
    sd_off = sd(mean_treg_off),
    mean_on = mean(mean_treg_on),
    mean_off = mean(mean_treg_off),
    mean = mean(mean_diff_on_off),
    min = min(mean_diff_on_off),
    max = max(mean_diff_on_off),
    variance = var(mean_diff_on_off),
    sd = sd(mean_diff_on_off),
    n_replicates = n()
  ) %>%
  arrange(desc(variance))

df_raw = df_raw %>% inner_join(var_df, by='param_set_id')
df_raw = df_raw %>% dplyr::mutate(high_var = ifelse(sd_on>1 & sd_off>1,1,0))
df_raw = df_raw %>% inner_join(df_params, by='param_set_id')


df_raw_low_var    = df_raw %>% dplyr::filter(high_var==0)
df_raw_high_var   = df_raw %>% dplyr::filter(high_var==1)

df_check_low_var = distinct(df_raw_low_var[c('param_set_id','mean','high_var')])
df_check_high_var = distinct(df_raw_high_var[c('param_set_id','mean','high_var')])

hist(df_check_low_var$mean) #!!!! -> THERE ARE TREGS BEING USEFUL IN THE LOW VARIANCE CASES! NEGATIVES ARE INSIGNIFICANT!
hist(df_check_high_var$mean) #!!!! -> THERE ARE TREGS BEING USEFUL IN THE LOW VARIANCE CASES! NEGATIVES ARE INSIGNIFICANT!

#### --- BE CAREFUL INTERPRETING HERE: THESE ARE THE "STABLE" CASES WITH LOW VARIATION
#### --- THEY MIGHT NOT AGREE WITH THE HISTOGRAMS FOR DECIDING LOW/HIGH VARIATION BECAUSE THAT'S A DIFFERENT THING!
#### --- 
df_check_low_var    = df_check_low_var %>% inner_join(df_params, by='param_set_id')
df_check_low_var_nz = df_check_low_var %>% dplyr::filter(mean>0)
# VERY STRONG CORRELATION -> WHY ARE TREGS MORE USEFUL WHEN activation_threshold_DAMPs is low? 
# Makes sense, that's when you need them?!
plot_param_vs_param(df_check_low_var_nz,'activation_threshold_DAMPs','mean') 
plot_param_vs_param(df_check_low_var_nz,'th_ROS_epith_recover','mean') # same goes for this too 

# ------------
plot_param_vs_param(df_raw, "sd_on", "mean_on")

# Filter for the two groups
df_plot = df_raw %>%
  filter(high_var %in% c(0, 1)) %>%
  mutate(high_var = factor(high_var, labels = c("Low Variance", "High Variance")))

library(ggplot2)

# vector of parameter names
param_names = c(
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
  "treg_discrimination_efficiency",
  "mean_on", # mean over reps
  "mean_off",# mean over reps
  "mean_treg_on",
  "mean_treg_off"
)

# loop over parameters
dir_name = './param_histograms'
dir.create(dir_name, showWarnings = FALSE)
for (param in param_names) {
  p = ggplot(df_plot, aes(x = .data[[param]], fill = high_var)) +
    geom_density(alpha = 0.5, bw = 0.1) +
    scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
    labs(x = param, title = paste("Density of", param)) +
    theme_minimal()
  ggsave(
    filename = paste0(dir_name, "/", param, "_density.png"),
    plot = p,
    width = 9,
    height = 6,
    dpi = 300
  )
}

df_plot_low_var   = df_plot %>% dplyr::filter(high_var=="Low Variance")
df_plot_high_var  = df_plot %>% dplyr::filter(high_var=="High Variance")

table(df_plot$high_var)

hist(df_plot_high_var$mean_treg_on)

sort(unique(round(df_plot_high_var$mean_treg_on,1)))
sort(unique(round(df_plot_low_var$mean_treg_on,1)))

saveRDS(df_plot,'df_plot_A15.rds')

