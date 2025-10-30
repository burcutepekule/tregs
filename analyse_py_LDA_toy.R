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

lda_data = data.frame(
  x_1 = runif(1000, min = -1, max = 1),
  x_2 = runif(1000, min = -1, max = 1)
)

# lda_data = lda_data %>% dplyr::mutate(class=ifelse(sign(x_1*x_2)<0 & abs(x_1*x_2)>1e-1,"N",
#                                                    ifelse(sign(x_1*x_2)>0 & abs(x_1*x_2)>1e-1,"P","0")))

lda_data = lda_data %>% dplyr::mutate(class=ifelse(sign(x_1+x_2)<0 & abs(x_1+x_2)>1e-1,"N",
                                                   ifelse(sign(x_1+x_2)>0 & abs(x_1+x_2)>1e-1,"P","0")))

table(lda_data$class)
                                      

# Fit LDA model
lda_model = lda(class ~ ., data = lda_data)
print(lda_model)

lda_pred = predict(lda_model)
lda_data$LD1 = lda_pred$x[,1]
lda_data$LD2 = lda_pred$x[,2]

lda_scores = as.data.frame(lda_pred$x)
lda_scores$class = lda_data$class

loadings = as.data.frame(lda_model$scaling)
loadings$feature = rownames(loadings)
loadings$LD1 = loadings$LD1 / max(abs(loadings$LD1)) * max(lda_scores$LD1)
loadings$LD2 = loadings$LD2 / max(abs(loadings$LD2)) * max(lda_scores$LD2)

# Compute vector length (Euclidean norm) of each feature's contribution
loadings$importance = sqrt(loadings$LD1^2 + loadings$LD2^2)

# Keep top N features (e.g., top 8)
top_features = loadings %>% filter(importance > 0.00*max(loadings$importance))

# Calculate percentages for each class_label
percentage_data = lda_scores %>%
  group_by(class) %>%
  summarise(
    count = n(),
    percentage = (n() / nrow(lda_scores)) * 100
  ) %>%
  mutate(
    label = paste0(class, "\n(", sprintf("%.2f", percentage), "%)")
  )

lda_scores_labeled = lda_scores %>% left_join(percentage_data, by = "class")

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
