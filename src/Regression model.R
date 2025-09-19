#### Dependencies ####
library(ggplot2)
library(dplyr)

# Load and clean up data, set global variables
source("src/import_data.R")   # defines colors_list
load("data/processed_data.RData")



#### Regression model ####

ggplot(data_combined, aes(x = Length, y = PC1)) +
  geom_point()


# Regression model with line #

ggplot(data_combined, aes(x = Length, y = PC1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) 


m <- lm(PC1 ~ Length, data = data_combined)

summary(m)


#### Correlation Analysis ####
# Overall Pearson correlation
cor_test <- cor.test(data_combined$PC1, data_combined$Length)
cat("Overall Pearson correlation:\n",
    "R =", round(cor_test$estimate, 3),
    "| p-value =", signif(cor_test$p.value, 3), "\n\n")


# ANOVA
anova_overall <- lm(PC1 ~ Length, data = data_combined)
summary(anova_overall)

# Plot
ggplot(data_combined, aes(x = Length, y = PC1, color = Region)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = colors_list$Region) +
  theme_minimal() +
  labs(
    title = "PC1 vs Length by Region",
    x = "Length [µm]",
    y = "PC1",
    color = "Region"
  )

