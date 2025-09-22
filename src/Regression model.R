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


confint(m)
coef(m)
ggplot(data_combined, aes(Length, PC1)) +
  geom_point() +
  geom_abline(aes(intercept = coef(m)[1], slope = coef(m)[2]),
              colour = "red")



#### ANOVA ####
n  <- aov(PC1 ~ Length, data = data_combined)
summary(n)



#### Pearson correlation ####
cor_test <- cor.test(data_combined$PC1, data_combined$Length)
cat("Overall Pearson correlation:\n",
    "R =", round(cor_test$estimate, 3),
    "| p-value =", signif(cor_test$p.value, 3), "\n\n")


