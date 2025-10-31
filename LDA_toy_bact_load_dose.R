rm(list=ls())
library(MASS)
library(ggplot2)
library(dplyr)

set.seed(123)

n         = 800
dose      = runif(n, 0 ,1)
bact_load = runif(n, 0, 1)
bact_res  = runif(n, 0, 1)

# works = ifelse(
#   (bact_load < 0.5 & dose > 0.5 & bact_res < 0.5) |
#     (bact_load < 0.5 & dose > 0.75 & (bact_res >= 0.5 & bact_res < 0.75)),
#   "works", "fails"
# )

works = ifelse(
  (bact_load < 0.5 & dose > 0.5 & bact_res < 0.5) |
    (bact_load < 0.5 & dose > 0.75 & bact_res >= 0.5 & bact_res < 0.75) |
    (bact_load > 0.7 & dose < 0.3),
  "works", "fails"
)

df_data = data.frame(dose, bact_load, bact_res, works = factor(works))

# Fit LDA
lda_model = lda(works ~ dose + bact_load + bact_res, data = df_data)
df_data$pred = predict(lda_model, df_data)$class
acc = mean(df_data$works == df_data$pred)
cat(sprintf("LDA accuracy: %.1f%%\n", 100*acc))

# Predict on grid for visualization
grid = expand.grid(
  dose = seq(0, 1, length = 200),
  bact_load = seq(0, 1, length = 200),
  bact_res = 0.5  # fix this dimension to visualize 2D slice
)

grid$pred = predict(lda_model, grid)$class

# Plot
ggplot() +
  geom_tile(data = grid, aes(x = dose, y = bact_load, fill = pred), alpha = 0.2) +
  geom_point(data = df_data, aes(x = dose, y = bact_load, color = works), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("fails" = "#d73027", "works" = "#1a9850")) +
  scale_fill_manual(values = c("fails" = "#fdae61", "works" = "#66bd63")) +
  labs(
    title = sprintf("LDA decision boundary (acc = %.1f%%)", 100*acc),
    x = "Dose", y = "Bacterial Load",
    fill = "LDA prediction", color = "True label"
  ) +
  theme_minimal()


# 1) Discriminant direction (LD1) in the 2D plane
w <- as.numeric(lda_model$scaling[c("dose","bact_load"), 1])
w <- w / sqrt(sum(w^2))              # unit length

center <- colMeans(df_data[, c("dose","bact_load")])
L <- 0.35                            # arrow half-length for display
seg_ld1 <- data.frame(
  x  = center[1] - L*w[1],
  y  = center[2] - L*w[2],
  xend = center[1] + L*w[1],
  yend = center[2] + L*w[2]
)

# 2) “Predictor vectors” via correlations with LD1 (structure matrix)
ld1_scores <- predict(lda_model)$x[,1]
cors <- cor(df_data[, c("dose","bact_load","bact_res")], ld1_scores)
cors2d <- cors[c("dose","bact_load")]               # only those we can draw in this plane
v <- as.numeric(cors2d / sqrt(sum(cors2d^2)))       # direction in the plot
L2 <- 0.28
seg_pred <- data.frame(
  var = c("dose","bact_load"),
  x = center[1], y = center[2],
  xend = center[1] + L2 * c(v[1], 0),
  yend = center[2] + L2 * c(0, v[2])
)

# Base decision map + points (your previous plot)
p <- ggplot() +
  geom_tile(data = grid, aes(x = dose, y = bact_load, fill = pred), alpha = 0.2) +
  geom_point(data = df_data, aes(x = dose, y = bact_load, color = works), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("fails" = "#d73027", "works" = "#1a9850")) +
  scale_fill_manual(values = c("fails" = "#fdae61", "works" = "#66bd63")) +
  labs(x = "Dose", y = "Bacterial Load", fill = "LDA prediction", color = "True label") +
  theme_minimal()

# Add LD1 direction and predictor arrows
p +
  geom_segment(data = seg_ld1, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.02, "npc")), linewidth = 1) +
  annotate("text", x = seg_ld1$xend, y = seg_ld1$yend, label = "LD1", vjust = -0.6, fontface = "bold") +
  geom_segment(data = seg_pred,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.018, "npc")), linetype = "dashed") +
  geom_text(data = seg_pred,
            aes(x = xend, y = yend, label = var), vjust = -0.6)
