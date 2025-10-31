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

set.seed(42)

# create simple dataset
n      = 500
gender = sample(c("man", "woman"), n, replace = TRUE)
age    = runif(n, 18, 60)

# conditional rule
# treatment works if woman <35 OR man >30
treatment_works = ifelse(
  (gender == "woman" & age < 35) | (gender == "man" & age > 30),
  "works", "fails"
)

df = data.frame(age, gender, treatment_works = factor(treatment_works))
table(df$treatment_works)

# encode gender numerically for LDA
df$gender_num = ifelse(df$gender == "man", 1, 0)

lda_model = lda(treatment_works ~ age + gender_num, data = df)
pred = predict(lda_model, df)$class
df$pred = pred

acc = mean(df$pred == df$treatment_works)
cat(sprintf("LDA accuracy: %.1f%%\n", 100 * acc))

grid = expand.grid(
  age = seq(min(age), max(age), length.out = 200),
  gender_num = seq(0, 1, length.out = 200)
)
grid$pred = predict(lda_model, grid)$class

ggplot() +
  geom_raster(data = grid, aes(age, gender_num, fill = pred), alpha = 0.3) +
  geom_point(data = df, aes(age, gender_num, color = treatment_works), size = 1.5) +
  scale_y_continuous(breaks = c(0, 1), labels = c("woman", "man")) +
  labs(
    title = "Conditional rule: treatment works if (woman <35) OR (man >30)",
    subtitle = "LDA draws a single linear boundary — misses conditional structure",
    x = "Age", y = "Gender"
  ) +
  theme_minimal()

