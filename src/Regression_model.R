#### Subsampling and Regression Analysis ####
library(ggplot2)
library(dplyr)
library(visreg)
library(tidyr)
library(purrr)

# Load processed data
set.seed(42)
source("src/import_data.R")
data_combined$Section <- as.factor(data_combined$Section)

#### Exploratory analysis ####

# Not sure if this needs to go into the manuscript

length_model <- lm(PC1 ~ Length, data = data_combined)
summary(length_model)
visreg(length_model, gg = TRUE)
length_Section_model <- lm(PC1 ~ Length * Section, data = data_combined)


# This saves the regression output so it can be exported to doc
Section_summary <- summary(length_Section_model)
save(Section_summary, file="data/regression_Section_summary.RData") #loaded in table_regression_Section.Rmd

# these models we do not show in the manuscript because they suffer from uneven # of observations
length_Section_model <- lm(PC1 ~ Length * Section, data = data_combined)
summary(length_Section_model)
visreg(length_Section_model, gg = TRUE, "Length", by="Section")
length_facies_model <- lm(PC1 ~ Length * Section + Section, data = data_combined)
summary(length_facies_model)
visreg(length_facies_model, "Length", by="Section", gg = TRUE, layout=c(3,1))

#### Subsampling ####

# Function to subsample n observations from each Section level
subsample_by_Section <- function(data, n_per_group = 20, seed = 42) {
  data %>%
    group_by(Section) %>%
    sample_n(size = min(n_per_group, n()), replace = FALSE) %>%
    ungroup()
}


# Run subsampling 
n_iterations <- 100
n_per_group <- 20
model_results <- list()
fitted_values_list <- list()
subsample_list <- list()

for (i in 1:n_iterations) {
  subsample_data <- subsample_by_Section(data_combined, n_per_group = n_per_group)
  subsample_list[[i]] <- subsample_data
  
  model <- lm(PC1 ~ Length * Section, data = subsample_data)
  
  model_results[[i]] <- list(
    iteration = i,
    model = model,
    coefficients = coef(model),
    r_squared = summary(model)$r.squared,
    adj_r_squared = summary(model)$adj.r.squared,
    sigma = summary(model)$sigma
  )
  
  fitted_data <- subsample_data %>%
    select(ID, Section, Length, PC1) %>%
    mutate(
      iteration = i,
      fitted = fitted(model),
      residual = residuals(model)
    )
  fitted_values_list[[i]] <- fitted_data
}

# Combine all fitted values
all_fitted <- bind_rows(fitted_values_list)

# Extract coefficients across iterations
all_coefficients <- do.call(rbind, lapply(model_results, function(x) x$coefficients))
coef_names <- colnames(all_coefficients)

# Coefficient statistics
coef_summary <- data.frame(
  coefficient = coef_names,
  mean = colMeans(all_coefficients, na.rm = TRUE),
  sd = apply(all_coefficients, 2, sd, na.rm = TRUE),
  lower_95 = apply(all_coefficients, 2, quantile, 0.025, na.rm = TRUE),
  upper_95 = apply(all_coefficients, 2, quantile, 0.975, na.rm = TRUE),
  lower_50 = apply(all_coefficients, 2, quantile, 0.25, na.rm = TRUE),
  upper_50 = apply(all_coefficients, 2, quantile, 0.75, na.rm = TRUE)
)

# Raw coefficients 
print(coef_summary, digits = 4)

# Calculate slopes for each Section in each iteration
slope_data <- data.frame()
Sections <- levels(data_combined$Section)

for (i in 1:n_iterations) {
  coefs <- model_results[[i]]$coefficients
  
  # Base slope (for reference Section - first level)
  base_slope <- coefs["Length"]
  
  for (sect in Sections) {
    if (sect == Sections[1]) {
      # Reference category gets base slope
      slope <- base_slope
    } else {
      # Other Sections get base + interaction
      interaction_name <- paste0("Length:Section", sect)
      if (interaction_name %in% names(coefs) && !is.na(coefs[interaction_name])) {
        slope <- base_slope + coefs[interaction_name]
      } else {
        slope <- base_slope
      }
    }
    slope_data <- rbind(slope_data, data.frame(
      iteration = i,
      Section = sect,
      slope = slope
    ))
  }
}

# Summarize slopes by Section
slope_summary <- slope_data %>%
  group_by(Section) %>%
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE),
    lower_95 = quantile(slope, 0.025, na.rm = TRUE),
    upper_95 = quantile(slope, 0.975, na.rm = TRUE),
    lower_50 = quantile(slope, 0.25, na.rm = TRUE),
    upper_50 = quantile(slope, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# Export for publication
save(slope_summary, file="data/slope_summary.RData") # used in slope_summary.Rmd

##### Line plot for individual iterations #####

p1 <- ggplot(all_fitted, aes(x = Length, y = fitted)) +
  geom_smooth(aes(group = interaction(Section, iteration), color = Section), 
              method = "lm", se = FALSE, alpha = 0.1, size = 0.3) +
  geom_smooth(aes(color = Section), method = "lm", se = TRUE, size = 1.2) +
  scale_color_manual(values = colors_list$Section) +
  labs(
    x = "Length [µm]",
    y = "PC1",
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )


p2 <- ggplot(all_fitted, aes(x = Length, y = fitted)) +
  geom_line(aes(group = iteration), alpha = 0.2, color = "gray50") +
  geom_smooth(method = "lm", color = "blue", fill = "blue", alpha = 0.3, size = 1) +
  geom_point(aes(y = PC1), alpha = 0.05, size = 0.5) +
  facet_wrap(~ Section, scales = "free_x", ncol = 3) +
  labs(
    x = "Length [µm]",
    y = "PC1",
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold")
  )

ggsave("supplementary_material/Fig.S4.jpg",p2, width=170, units="mm", dpi = 300)

##### Visualize slope distributions with confidence intervals #####

p3 <- ggplot(slope_data, aes(x = Section, y = slope, fill = Section)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_point(data = slope_summary, aes(y = mean_slope), 
             size = 3, color = "black", shape = 18) +
  geom_errorbar(data = slope_summary, 
                aes(y = mean_slope, ymin = lower_95, ymax = upper_95),
                width = 0.1, color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = colors_list$Section) +
  labs(
    x = "Section",
    y = "Slope (Length effect on PC1)",
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
  )
ggsave("figs/Fig.9.jpg",p3, width=170, units="mm", dpi = 300)

##### Compare slope confidence intervals #####

# a different version of the plot, not used in the end
p4 <- ggplot(slope_summary, aes(x = Section, y = mean_slope, color = Section)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, size = 1) +
  geom_errorbar(aes(ymin = lower_50, ymax = upper_50), width = 0.1, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = colors_list$Section) +
  labs(
    x = "Section",
    y = "Slope (Length effect on PC1)",
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
  )

# Thick bars: 50% CI, Thin bars: 95% CI
p4

##### R² #####

# R-squared distribution
r_squared_values <- sapply(model_results, function(x) x$r_squared)
# Mean R²
round(mean(r_squared_values), 4)
# SD R²
round(sd(r_squared_values), 4)
# 95% CI R²
round(quantile(r_squared_values, 0.025), 4) 
round(quantile(r_squared_values, 0.975), 4)

p5 <- ggplot(data.frame(r_squared = r_squared_values), aes(x = r_squared)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = mean(r_squared_values), color = "red", 
             linetype = "dashed", size = 1) +
  geom_vline(xintercept = quantile(r_squared_values, c(0.025, 0.975)), 
             color = "red", linetype = "dotted", size = 0.5) +
  labs(
    x = "R²",
    y = "Frequency",
  ) +
  theme_minimal() 

ggsave("supplementary_material/Fig.S6.jpg",p5, width=170, units="mm", dpi = 300) # Dashed line: mean, Dotted lines: 95% CI

##### ANOVA for interactions #####

# Variance components: how much variation is due to Section vs length?
for (i in 1:n_iterations) {
  anova_result <- anova(model_results[[i]]$model)
}

# Test if Section interactions are consistently significant
interaction_pvalues <- sapply(model_results, function(x) {
  anova_result <- anova(x$model)
  # p-value for Length:Section interaction
  if ("Length:Section" %in% rownames(anova_result)) {
    anova_result["Length:Section", "Pr(>F)"]
  } else {
    NA
  }
})

# Proportion of iterations with p < 0.05:
mean(interaction_pvalues < 0.05)/100
# Median p-value:
median(interaction_pvalues)

