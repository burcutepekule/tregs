rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(mgcv)
library(factoextra)
library(cluster)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")
source("~/Desktop/tregs/A013_datazanalyse_sterile_1_rnd_0.R")
df_plot = readRDS('df_plot_A13.rds') # comes from 
# # ---------- Conditional Inference Trees (CTree): These use statistical tests for splitting and handle interactions better:
library(partykit)

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
  "treg_discrimination_efficiency"
)

df_clustering = distinct(df_plot[c('high_var',param_names)])

df_clustering = df_clustering %>%
  mutate(dominant_outcome = factor(high_var, levels = c("High Variance","Low Variance")))


library(rpart)
library(rpart.plot)

# Calculate class weights (inverse of frequency)
class_weights = 1 / prop.table(table(df_clustering$dominant_outcome))

cart_purity = rpart(
  dominant_outcome ~ .,
  data = df_clustering %>% select(dominant_outcome, all_of(param_names)),
  method = "class",
  weights = ifelse(df_clustering$dominant_outcome == "Low Variance", 
                   class_weights["Low Variance"],
                   class_weights["High Variance"]),
  control = rpart.control(
    cp = 0.02,
    minsplit = 10,
    minbucket = 5,
    maxdepth = 10
  )
)

# Save a very high resolution plot to see all details
png("tree_high_res.png", width = 7000, height = 3000, res = 300)
# rpart.plot(cart_purity, 
#            type = 4, 
#            extra = 104,
#            under = FALSE,
#            faclen = 0,
#            cex = 0.6,  # Smaller text to fit more
#            tweak = 1.2)
# 

rpart.plot(
  cart_purity,
  type = 4,
  extra = 104,
  under = TRUE,      # put labels underneath
  faclen = 0,
  cex = 0.6,
  tweak = 1.2,
  fallen.leaves = TRUE,
  nn = TRUE          # <— shows node numbers on the plot
)
dev.off()


# Get full frame info
frame_info <- cart_purity$frame
# Add node IDs
frame_info$node_id <- as.numeric(row.names(frame_info))

terminal_nodes = frame_info %>%
  dplyr::filter(var == "<leaf>") %>%
  dplyr::mutate(
    Node_ID = as.numeric(row.names(.)),
    Predicted_Class = ifelse(yval == 1, "High Variance", "Low Variance"),
    N_Samples = n,
    P_HighVar = yval2[, 4],
    P_LowVar  = yval2[, 5],
    Node_Percent = yval2[, 6]
  ) %>%
  dplyr::select(Node_ID, Predicted_Class, N_Samples, P_HighVar, P_LowVar, Node_Percent)


