rm(list=ls())
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(zoo)
library(slider)
library(purrr)
library(tibble)
library(RcppRoll)
library(MASS)
library(ggrepel)
library(randomForest)
library(vip)
library(pdp)
df_short_merged_with_params = readRDS('/Users/burcutepekule/Desktop/tregs/df_short_merged_with_params.rds')

param_cols = colnames(df_short_merged_with_params)[8:31]
# Prepare data
lda_data = df_short_merged_with_params %>%
  dplyr::select(class_code, all_of(param_cols)) %>%
  filter(!is.na(class_code))  # Remove any NA classes

features = lda_data %>% dplyr::select(-class_code)
features = as.data.frame(scale(features))

# Fit LDA model
lda_model = lda(class_code ~ ., data = cbind(class_code = lda_data$class_code, features))
print(lda_model)

lda_pred = predict(lda_model)
lda_data$LD1 = lda_pred$x[,1]
lda_data$LD2 = lda_pred$x[,2]

lda_scores = as.data.frame(lda_pred$x)
lda_scores$class_code = lda_data$class_code

loadings = as.data.frame(lda_model$scaling)
loadings$feature = rownames(loadings)
loadings$LD1 = loadings$LD1 / max(abs(loadings$LD1)) * max(lda_scores$LD1)
loadings$LD2 = loadings$LD2 / max(abs(loadings$LD2)) * max(lda_scores$LD2)

# Compute vector length (Euclidean norm) of each feature's contribution
loadings$importance = sqrt(loadings$LD1^2 + loadings$LD2^2)

# Keep top N features (e.g., top 8)
top_features = loadings %>% filter(importance > 0.00*max(loadings$importance))

lda_scores = lda_scores %>% mutate(class_label = case_when(
  class_code == "0000" ~ "Drift", #"No effect",
  class_code == "1000" ~ "Favorable", #Sterile better only",
  class_code == "0100" ~ "Favorable", #"Pathogenic better only",
  class_code == "0010" ~ "Unfavorable", #"Sterile worse only",
  class_code == "0001" ~ "Unfavorable", #"Pathogenic worse only",
  class_code == "1100" ~ "Favorable", #"Both better",
  class_code == "0011" ~ "Unfavorable", #"Both worse",
  class_code == "1001" ~ "Unfavorable", # "Sterile better, pathogenic worse",
  class_code == "0110" ~ "Unfavorable", # "Pathogenic better, sterile worse",
  TRUE ~ "Other"
))

# Calculate percentages for each class_label
percentage_data = lda_scores %>%
  group_by(class_label) %>%
  summarise(
    count = n(),
    percentage = (n() / nrow(lda_scores)) * 100
  ) %>%
  mutate(
    label = paste0(class_label, "\n(", sprintf("%.2f", percentage), "%)")
  )

lda_scores_labeled = lda_scores %>% left_join(percentage_data, by = "class_label")

library(RColorBrewer)
custom_colors = setNames(
  c("#AFAAB9","blue","red"),
  percentage_data$label
)

top_features       = top_features[order(-top_features$importance), ]
top_features$index = seq_len(nrow(top_features))


p_lda = ggplot(lda_scores_labeled, aes(x = LD1, y = LD2, color = label)) +
  geom_point(alpha = 0.4, size = 2) +
  stat_ellipse(type = "norm", level = 0.95, size = 1, linetype = "dashed") +
  geom_segment(data = top_features,
               aes(x = 0, y = 0, xend = LD1, yend = LD2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 1.1) +
  geom_text_repel(data = top_features,
                  aes(x = LD1, y = LD2, label = index),
                  color = "black", size = 4.5, fontface = "bold", 
                  force = 5, max.overlaps = Inf, segment.color = NA) +
  scale_color_manual(values = custom_colors) +
  labs(title = "LDA projection (2D) with top features, pathogenic injury",
       x = "LD1", y = "LD2", color = "Treg effect class") +
  theme_minimal(base_size = 14)

feature_table = tableGrob(
  data.frame(Index = top_features$index, Feature = top_features$feature, LD1 = round(top_features$LD1,2), LD2=round(top_features$LD2,2)),
  rows = NULL, theme = ttheme_minimal(base_size = 10)
)

final_plot = ggdraw() +
  draw_plot(p_lda, 0, 0, 0.75, 1) +
  draw_grob(feature_table, 0.76, 0.1, 0.24, 0.8)  # adjust x, y, width, height as needed

print(final_plot)
