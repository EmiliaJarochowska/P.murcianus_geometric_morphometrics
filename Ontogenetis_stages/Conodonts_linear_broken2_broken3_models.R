# ============================================================
# LINEAR REGRESSION + BROKEN-STICK MODELS WITH 1 AND 2 BREAKPOINTS
# ============================================================
#
# Models:
#   1) Linear: one straight line
#   2) Broken-stick with 1 breakpoint: two connected lines
#   3) Broken-stick with 2 breakpoints: three connected lines
#
# The main comparison is based on AICc. The script also reports
# AIC, BIC, RMSE, R-squared, Akaike weights, equations, slopes,
# breakpoints, and publication-ready figures.
# ============================================================


# 1. Packages -------------------------------------------------

packages <- c(
  "readxl",
  "ggplot2",
  "segmented",
  "patchwork"
)

missing_packages <- setdiff(
  packages,
  rownames(installed.packages())
)

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(
  lapply(
    packages,
    library,
    character.only = TRUE
  )
)

set.seed(123)


# 2. Input file ------------------------------------------------

file_path <- paste0(
  "C:/Users/Humbe/OneDrive - University of Bristol/",
  "Desktop/Conodontos_Carlos/data_combined.xlsx"
)

output_dir <- dirname(file_path)


# 3. Import and prepare the data -------------------------------

raw_data <- readxl::read_excel(
  path = file_path,
  sheet = 1
)

if (ncol(raw_data) < 3) {
  stop(
    "The Excel file must contain at least three columns: ",
    "ID, Length, and Denticles."
  )
}

dat <- as.data.frame(
  raw_data[, 1:3]
)

names(dat) <- c(
  "ID",
  "Length",
  "Denticles"
)

dat$Length <- suppressWarnings(
  as.numeric(dat$Length)
)

dat$Denticles <- suppressWarnings(
  as.numeric(dat$Denticles)
)

dat <- dat[
  complete.cases(
    dat[, c("Length", "Denticles")]
  ),
]

if (nrow(dat) < 15) {
  warning(
    "There are few observations. Models with two breakpoints ",
    "may be unstable or overfitted."
  )
}

if (length(unique(dat$Length)) < 8) {
  stop(
    "There are not enough distinct Length values to fit ",
    "a three-segment model."
  )
}

cat(
  "Number of observations:",
  nrow(dat),
  "\n"
)

cat(
  "Length range:",
  min(dat$Length),
  "-",
  max(dat$Length),
  "\n\n"
)


# 4. Simple linear regression ---------------------------------

model_linear <- lm(
  Denticles ~ Length,
  data = dat
)

cat(
  "SIMPLE LINEAR REGRESSION\n"
)

print(
  summary(model_linear)
)


# 5. Helper functions ------------------------------------------

# Count observations in the segments defined by a set of breakpoints.
segment_counts <- function(
  x,
  breakpoints
) {
  breakpoints <- sort(
    as.numeric(breakpoints)
  )

  intervals <- cut(
    x,
    breaks = c(
      -Inf,
      breakpoints,
      Inf
    ),
    right = TRUE
  )

  as.integer(
    table(intervals)
  )
}


# Check that a segmented model is usable and that each segment
# contains a minimum number of observations.
is_valid_segmented_model <- function(
  model,
  expected_breakpoints,
  x,
  min_n_per_segment
) {
  if (inherits(model, "try-error")) {
    return(FALSE)
  }

  if (
    is.null(model$psi) ||
    nrow(model$psi) != expected_breakpoints
  ) {
    return(FALSE)
  }

  breakpoints <- sort(
    as.numeric(
      model$psi[, "Est."]
    )
  )

  if (
    any(!is.finite(breakpoints)) ||
    any(diff(breakpoints) <= 0)
  ) {
    return(FALSE)
  }

  counts <- segment_counts(
    x,
    breakpoints
  )

  length(counts) == expected_breakpoints + 1 &&
    all(counts >= min_n_per_segment)
}


# Fit many starting solutions and retain the valid model
# with the lowest AIC.
fit_segmented_multistart <- function(
  base_model,
  starting_values,
  expected_breakpoints,
  x,
  min_n_per_segment,
  n_boot = 10
) {
  model_attempts <- lapply(
    starting_values,
    function(initial_psi) {
      try(
        segmented::segmented(
          obj = base_model,
          seg.Z = ~ Length,
          psi = list(
            Length = initial_psi
          ),
          control = segmented::seg.control(
            display = FALSE,
            it.max = 300,
            n.boot = n_boot,
            seed = 123
          )
        ),
        silent = TRUE
      )
    }
  )

  valid_models <- Filter(
    function(model) {
      is_valid_segmented_model(
        model = model,
        expected_breakpoints = expected_breakpoints,
        x = x,
        min_n_per_segment = min_n_per_segment
      )
    },
    model_attempts
  )

  if (length(valid_models) == 0) {
    return(NULL)
  }

  model_AICs <- vapply(
    valid_models,
    AIC,
    numeric(1)
  )

  valid_models[[which.min(model_AICs)]]
}


# Extract the equations, slopes, intersections, and approximate
# Wald intervals from a segmented model.
extract_segmented_parameters <- function(
  model,
  variable_name = "Length"
) {
  psi_table <- as.data.frame(
    model$psi
  )

  psi_table <- psi_table[
    order(psi_table$Est.),
    ,
    drop = FALSE
  ]

  breakpoints <- as.numeric(
    psi_table$Est.
  )

  breakpoint_se <- as.numeric(
    psi_table$St.Err
  )

  breakpoint_ci_lower <-
    breakpoints +
    qnorm(0.025) * breakpoint_se

  breakpoint_ci_upper <-
    breakpoints +
    qnorm(0.975) * breakpoint_se

  model_coef <- coef(model)

  base_intercept <- unname(
    model_coef["(Intercept)"]
  )

  base_slope <- unname(
    model_coef[variable_name]
  )

  hinge_names <- grep(
    pattern = paste0(
      "^U[0-9]+\\.",
      variable_name,
      "$"
    ),
    x = names(model_coef),
    value = TRUE
  )

  hinge_numbers <- as.integer(
    sub(
      paste0(
        "^U([0-9]+)\\.",
        variable_name,
        "$"
      ),
      "\\1",
      hinge_names
    )
  )

  hinge_names <- hinge_names[
    order(hinge_numbers)
  ]

  slope_changes <- unname(
    model_coef[hinge_names]
  )

  if (
    length(slope_changes) !=
    length(breakpoints)
  ) {
    stop(
      "The number of slope-change coefficients does not match ",
      "the number of estimated breakpoints."
    )
  }

  slopes <- c(
    base_slope,
    base_slope + cumsum(
      slope_changes
    )
  )

  intercepts <- numeric(
    length(slopes)
  )

  intercepts[1] <- base_intercept

  if (length(breakpoints) >= 1) {
    for (
      segment_index in
      2:length(slopes)
    ) {
      intercepts[segment_index] <-
        base_intercept -
        sum(
          slope_changes[
            seq_len(
              segment_index - 1
            )
          ] *
            breakpoints[
              seq_len(
                segment_index - 1
              )
            ]
        )
    }
  }

  intersection_y <- as.numeric(
    predict(
      model,
      newdata = data.frame(
        Length = breakpoints
      )
    )
  )

  list(
    breakpoints = breakpoints,
    breakpoint_se = breakpoint_se,
    breakpoint_ci_lower = breakpoint_ci_lower,
    breakpoint_ci_upper = breakpoint_ci_upper,
    slope_changes = slope_changes,
    slopes = slopes,
    intercepts = intercepts,
    intersection_y = intersection_y
  )
}


# Calculate comparable statistics for all candidate models.
get_model_statistics <- function(
  model,
  model_name,
  response
) {
  model_logLik <- logLik(model)

  k <- attr(
    model_logLik,
    "df"
  )

  n <- nobs(model)

  model_AIC <- AIC(model)

  model_AICc <- if (
    n > k + 1
  ) {
    model_AIC +
      (
        2 * k * (k + 1)
      ) /
      (
        n - k - 1
      )
  } else {
    NA_real_
  }

  residual_values <- residuals(model)

  rss <- sum(
    residual_values^2
  )

  rmse <- sqrt(
    mean(
      residual_values^2
    )
  )

  total_ss <- sum(
    (
      response -
      mean(response)
    )^2
  )

  r_squared <- 1 -
    rss / total_ss

  data.frame(
    Model = model_name,
    N = n,
    K = k,
    logLik = as.numeric(
      model_logLik
    ),
    AIC = model_AIC,
    AICc = model_AICc,
    BIC = BIC(model),
    RMSE = rmse,
    R_squared = r_squared
  )
}


# 6. Broken-stick with 1 breakpoint / 2 lines -----------------

one_break_starting_values <- lapply(
  unique(
    as.numeric(
      quantile(
        dat$Length,
        probs = c(
          0.20,
          0.30,
          0.40,
          0.50,
          0.60,
          0.70,
          0.80
        ),
        na.rm = TRUE
      )
    )
  ),
  function(x) x
)

min_n_one_break <- max(
  10,
  ceiling(
    0.05 * nrow(dat)
  )
)

model_broken_2lines <- fit_segmented_multistart(
  base_model = model_linear,
  starting_values = one_break_starting_values,
  expected_breakpoints = 1,
  x = dat$Length,
  min_n_per_segment = min_n_one_break,
  n_boot = 10
)

if (is.null(model_broken_2lines)) {
  stop(
    "The broken-stick model with one breakpoint did not ",
    "produce a valid solution."
  )
}

parameters_2lines <- extract_segmented_parameters(
  model_broken_2lines
)

cat(
  "\n\nBROKEN-STICK: 1 BREAKPOINT / 2 LINES\n"
)

print(
  summary(model_broken_2lines)
)


# 7. Broken-stick with 2 breakpoints / 3 lines ----------------

# Candidate starting pairs are generated from internal quantiles.
# Pairs that are too close together are removed.
start_quantile_probabilities <- seq(
  0.15,
  0.85,
  by = 0.10
)

start_quantiles <- unique(
  as.numeric(
    quantile(
      dat$Length,
      probs = start_quantile_probabilities,
      na.rm = TRUE
    )
  )
)

starting_pair_matrix <- t(
  combn(
    start_quantiles,
    2
  )
)

minimum_starting_separation <-
  0.12 *
  diff(
    range(dat$Length)
  )

starting_pair_matrix <- starting_pair_matrix[
  (
    starting_pair_matrix[, 2] -
    starting_pair_matrix[, 1]
  ) >= minimum_starting_separation,
  ,
  drop = FALSE
]

two_break_starting_values <- lapply(
  seq_len(
    nrow(starting_pair_matrix)
  ),
  function(i) {
    as.numeric(
      starting_pair_matrix[i, ]
    )
  }
)

min_n_two_breaks <- max(
  15,
  ceiling(
    0.05 * nrow(dat)
  )
)

model_broken_3lines <- fit_segmented_multistart(
  base_model = model_linear,
  starting_values = two_break_starting_values,
  expected_breakpoints = 2,
  x = dat$Length,
  min_n_per_segment = min_n_two_breaks,
  n_boot = 10
)

if (is.null(model_broken_3lines)) {
  stop(
    "The broken-stick model with two breakpoints did not ",
    "produce a valid solution. Try reducing min_n_two_breaks ",
    "or inspect the data distribution."
  )
}

parameters_3lines <- extract_segmented_parameters(
  model_broken_3lines
)

cat(
  "\n\nBROKEN-STICK: 2 BREAKPOINTS / 3 LINES\n"
)

print(
  summary(model_broken_3lines)
)


# 8. Equations and breakpoint tables ---------------------------

linear_coef <- coef(
  model_linear
)

linear_intercept <- unname(
  linear_coef["(Intercept)"]
)

linear_slope <- unname(
  linear_coef["Length"]
)


parameter_table_2lines <- data.frame(
  Model = "Broken-stick: 2 lines",
  Segment = c(
    "Line 1",
    "Line 2"
  ),
  Intercept = parameters_2lines$intercepts,
  Slope = parameters_2lines$slopes
)

breakpoint_table_2lines <- data.frame(
  Model = "Broken-stick: 2 lines",
  Breakpoint = 1,
  Length = parameters_2lines$breakpoints,
  Denticles = parameters_2lines$intersection_y,
  Standard_error = parameters_2lines$breakpoint_se,
  CI_95_lower = parameters_2lines$breakpoint_ci_lower,
  CI_95_upper = parameters_2lines$breakpoint_ci_upper
)


parameter_table_3lines <- data.frame(
  Model = "Broken-stick: 3 lines",
  Segment = c(
    "Line 1",
    "Line 2",
    "Line 3"
  ),
  Intercept = parameters_3lines$intercepts,
  Slope = parameters_3lines$slopes
)

breakpoint_table_3lines <- data.frame(
  Model = "Broken-stick: 3 lines",
  Breakpoint = seq_along(
    parameters_3lines$breakpoints
  ),
  Length = parameters_3lines$breakpoints,
  Denticles = parameters_3lines$intersection_y,
  Standard_error = parameters_3lines$breakpoint_se,
  CI_95_lower = parameters_3lines$breakpoint_ci_lower,
  CI_95_upper = parameters_3lines$breakpoint_ci_upper
)


cat(
  "\n\nLINEAR MODEL EQUATION\n"
)

cat(
  sprintf(
    "Denticles = %.6f %+.6f x Length\n",
    linear_intercept,
    linear_slope
  )
)


cat(
  "\nTWO-LINE BROKEN-STICK EQUATIONS\n"
)

for (
  i in seq_len(
    nrow(parameter_table_2lines)
  )
) {
  cat(
    sprintf(
      "%s: Denticles = %.6f %+.6f x Length\n",
      parameter_table_2lines$Segment[i],
      parameter_table_2lines$Intercept[i],
      parameter_table_2lines$Slope[i]
    )
  )
}

print(
  breakpoint_table_2lines,
  row.names = FALSE
)


cat(
  "\nTHREE-LINE BROKEN-STICK EQUATIONS\n"
)

for (
  i in seq_len(
    nrow(parameter_table_3lines)
  )
) {
  cat(
    sprintf(
      "%s: Denticles = %.6f %+.6f x Length\n",
      parameter_table_3lines$Segment[i],
      parameter_table_3lines$Intercept[i],
      parameter_table_3lines$Slope[i]
    )
  )
}

print(
  breakpoint_table_3lines,
  row.names = FALSE
)


# 9. Compare the three models ---------------------------------

model_comparison <- rbind(
  get_model_statistics(
    model_linear,
    "Linear",
    dat$Denticles
  ),
  get_model_statistics(
    model_broken_2lines,
    "Broken-stick: 2 lines",
    dat$Denticles
  ),
  get_model_statistics(
    model_broken_3lines,
    "Broken-stick: 3 lines",
    dat$Denticles
  )
)

model_comparison$Delta_AIC <-
  model_comparison$AIC -
  min(
    model_comparison$AIC
  )

model_comparison$AIC_weight <-
  exp(
    -0.5 *
    model_comparison$Delta_AIC
  ) /
  sum(
    exp(
      -0.5 *
      model_comparison$Delta_AIC
    )
  )

model_comparison$Delta_AICc <-
  model_comparison$AICc -
  min(
    model_comparison$AICc,
    na.rm = TRUE
  )

model_comparison$AICc_weight <-
  exp(
    -0.5 *
    model_comparison$Delta_AICc
  ) /
  sum(
    exp(
      -0.5 *
      model_comparison$Delta_AICc
    )
  )

model_comparison$Delta_BIC <-
  model_comparison$BIC -
  min(
    model_comparison$BIC
  )

model_comparison <- model_comparison[
  order(
    model_comparison$AICc
  ),
]

rownames(
  model_comparison
) <- NULL

cat(
  "\n\nMODEL COMPARISON\n"
)

print(
  model_comparison,
  digits = 5,
  row.names = FALSE
)

best_model <- model_comparison$Model[1]

cat(
  "\nBest-supported model according to AICc:",
  best_model,
  "\n"
)


# Specific comparison between the two segmented alternatives.
aicc_2lines <- model_comparison$AICc[
  model_comparison$Model ==
  "Broken-stick: 2 lines"
]

aicc_3lines <- model_comparison$AICc[
  model_comparison$Model ==
  "Broken-stick: 3 lines"
]

delta_AICc_3_vs_2 <-
  aicc_3lines -
  aicc_2lines

cat(
  sprintf(
    paste0(
      "\nAICc(3 lines) - AICc(2 lines) = %.3f\n",
      "Negative values favour three lines; ",
      "positive values favour two lines.\n"
    ),
    delta_AICc_3_vs_2
  )
)

if (delta_AICc_3_vs_2 <= -2) {
  cat(
    "Interpretation: the third line receives meaningful AICc support.\n"
  )
} else if (delta_AICc_3_vs_2 >= 2) {
  cat(
    "Interpretation: the additional breakpoint is not justified; ",
    "prefer the simpler two-line broken-stick.\n"
  )
} else {
  cat(
    "Interpretation: the two segmented models have similar AICc support; ",
    "prefer the simpler model unless the third phase has a strong ",
    "independent biological justification.\n"
  )
}


# 10. Predictions ----------------------------------------------

prediction_grid <- data.frame(
  Length = seq(
    from = min(dat$Length),
    to = max(dat$Length),
    length.out = 600
  )
)

prediction_grid$Linear <- predict(
  model_linear,
  newdata = prediction_grid
)

prediction_grid$Broken_stick_2lines <- predict(
  model_broken_2lines,
  newdata = prediction_grid
)

prediction_grid$Broken_stick_3lines <- predict(
  model_broken_3lines,
  newdata = prediction_grid
)


# 11. Individual model plots ----------------------------------

plot_linear <- ggplot(
  dat,
  aes(
    x = Length,
    y = Denticles
  )
) +
  geom_point(
    size = 2.2,
    alpha = 0.55
  ) +
  geom_line(
    data = prediction_grid,
    aes(
      x = Length,
      y = Linear
    ),
    linewidth = 1
  ) +
  labs(
    title = "Linear regression",
    subtitle = sprintf(
      "AICc = %.2f; R² = %.3f",
      model_comparison$AICc[
        model_comparison$Model ==
        "Linear"
      ],
      model_comparison$R_squared[
        model_comparison$Model ==
        "Linear"
      ]
    ),
    x = "Element length",
    y = "Number of denticles on the blade"
  ) +
  theme_classic(
    base_size = 12
  )


plot_broken_2lines <- ggplot(
  dat,
  aes(
    x = Length,
    y = Denticles
  )
) +
  geom_point(
    size = 2.2,
    alpha = 0.55
  ) +
  geom_line(
    data = prediction_grid,
    aes(
      x = Length,
      y = Broken_stick_2lines
    ),
    linewidth = 1
  ) +
  geom_vline(
    xintercept =
      parameters_2lines$breakpoints,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    data = breakpoint_table_2lines,
    aes(
      x = Length,
      y = Denticles
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  labs(
    title = "Broken-stick: two lines",
    subtitle = sprintf(
      "Breakpoint = %.2f; AICc = %.2f; R² = %.3f",
      parameters_2lines$breakpoints,
      model_comparison$AICc[
        model_comparison$Model ==
        "Broken-stick: 2 lines"
      ],
      model_comparison$R_squared[
        model_comparison$Model ==
        "Broken-stick: 2 lines"
      ]
    ),
    x = "Element length",
    y = "Number of denticles on the blade"
  ) +
  theme_classic(
    base_size = 12
  )


plot_broken_3lines <- ggplot(
  dat,
  aes(
    x = Length,
    y = Denticles
  )
) +
  geom_point(
    size = 2.2,
    alpha = 0.55
  ) +
  geom_line(
    data = prediction_grid,
    aes(
      x = Length,
      y = Broken_stick_3lines
    ),
    linewidth = 1
  ) +
  geom_vline(
    xintercept =
      parameters_3lines$breakpoints,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    data = breakpoint_table_3lines,
    aes(
      x = Length,
      y = Denticles
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  labs(
    title = "Broken-stick: three lines",
    subtitle = sprintf(
      paste0(
        "Breakpoints = %.2f and %.2f; ",
        "AICc = %.2f; R² = %.3f"
      ),
      parameters_3lines$breakpoints[1],
      parameters_3lines$breakpoints[2],
      model_comparison$AICc[
        model_comparison$Model ==
        "Broken-stick: 3 lines"
      ],
      model_comparison$R_squared[
        model_comparison$Model ==
        "Broken-stick: 3 lines"
      ]
    ),
    x = "Element length",
    y = "Number of denticles on the blade"
  ) +
  theme_classic(
    base_size = 12
  )


combined_plot <- (
  plot_linear |
  plot_broken_2lines |
  plot_broken_3lines
) +
  patchwork::plot_annotation(
    title = paste(
      "Ontogenetic relationship between element length",
      "and denticle number"
    )
  )

print(
  combined_plot
)


# 12. Overlay of all three models ------------------------------

comparison_predictions <- rbind(
  data.frame(
    Length = prediction_grid$Length,
    Prediction = prediction_grid$Linear,
    Model = "Linear"
  ),
  data.frame(
    Length = prediction_grid$Length,
    Prediction =
      prediction_grid$Broken_stick_2lines,
    Model = "Broken-stick: 2 lines"
  ),
  data.frame(
    Length = prediction_grid$Length,
    Prediction =
      prediction_grid$Broken_stick_3lines,
    Model = "Broken-stick: 3 lines"
  )
)

comparison_predictions$Model <- factor(
  comparison_predictions$Model,
  levels = c(
    "Linear",
    "Broken-stick: 2 lines",
    "Broken-stick: 3 lines"
  )
)

overlay_plot <- ggplot(
  dat,
  aes(
    x = Length,
    y = Denticles
  )
) +
  geom_point(
    size = 2.2,
    alpha = 0.50
  ) +
  geom_line(
    data = comparison_predictions,
    aes(
      x = Length,
      y = Prediction,
      linetype = Model
    ),
    linewidth = 1.05
  ) +
  scale_linetype_manual(
    values = c(
      "Linear" = "dashed",
      "Broken-stick: 2 lines" = "solid",
      "Broken-stick: 3 lines" = "dotdash"
    )
  ) +
  labs(
    title = "Comparison of the three candidate models",
    subtitle = paste0(
      "Best-supported model according to AICc: ",
      best_model
    ),
    x = "Element length",
    y = "Number of denticles on the blade",
    linetype = "Model"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "top"
  )

print(
  overlay_plot
)


# 13. AICc comparison plot ------------------------------------

aicc_plot_data <- model_comparison

aicc_plot_data$Model <- factor(
  aicc_plot_data$Model,
  levels = rev(
    aicc_plot_data$Model
  )
)

aicc_plot <- ggplot(
  aicc_plot_data,
  aes(
    x = Model,
    y = Delta_AICc
  )
) +
  geom_col(
    width = 0.65
  ) +
  coord_flip() +
  geom_hline(
    yintercept = 2,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    title = "Relative support for the candidate models",
    subtitle = "Lower Delta AICc indicates greater support",
    x = NULL,
    y = expression(Delta * AIC[c])
  ) +
  theme_classic(
    base_size = 12
  )

print(
  aicc_plot
)


# 14. Save figures and tables ---------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "linear_broken2_broken3_panels.png"
  ),
  plot = combined_plot,
  width = 16,
  height = 5.5,
  dpi = 400
)

ggsave(
  filename = file.path(
    output_dir,
    "linear_broken2_broken3_panels.pdf"
  ),
  plot = combined_plot,
  width = 16,
  height = 5.5
)

ggsave(
  filename = file.path(
    output_dir,
    "linear_broken2_broken3_overlay.png"
  ),
  plot = overlay_plot,
  width = 8,
  height = 5.8,
  dpi = 400
)

ggsave(
  filename = file.path(
    output_dir,
    "linear_broken2_broken3_overlay.pdf"
  ),
  plot = overlay_plot,
  width = 8,
  height = 5.8
)

ggsave(
  filename = file.path(
    output_dir,
    "model_comparison_delta_AICc.png"
  ),
  plot = aicc_plot,
  width = 7,
  height = 4.5,
  dpi = 400
)

ggsave(
  filename = file.path(
    output_dir,
    "model_comparison_delta_AICc.pdf"
  ),
  plot = aicc_plot,
  width = 7,
  height = 4.5
)


write.csv(
  model_comparison,
  file = file.path(
    output_dir,
    "model_comparison_3_models.csv"
  ),
  row.names = FALSE
)

write.csv(
  parameter_table_2lines,
  file = file.path(
    output_dir,
    "broken_stick_2lines_equations.csv"
  ),
  row.names = FALSE
)

write.csv(
  breakpoint_table_2lines,
  file = file.path(
    output_dir,
    "broken_stick_2lines_breakpoint.csv"
  ),
  row.names = FALSE
)

write.csv(
  parameter_table_3lines,
  file = file.path(
    output_dir,
    "broken_stick_3lines_equations.csv"
  ),
  row.names = FALSE
)

write.csv(
  breakpoint_table_3lines,
  file = file.path(
    output_dir,
    "broken_stick_3lines_breakpoints.csv"
  ),
  row.names = FALSE
)

write.csv(
  prediction_grid,
  file = file.path(
    output_dir,
    "model_predictions_3_models.csv"
  ),
  row.names = FALSE
)


capture.output(
  {
    cat(
      "SIMPLE LINEAR REGRESSION\n\n"
    )

    print(
      summary(model_linear)
    )

    cat(
      "\n\nBROKEN-STICK: 1 BREAKPOINT / 2 LINES\n\n"
    )

    print(
      summary(model_broken_2lines)
    )

    cat(
      "\n\nBROKEN-STICK: 2 BREAKPOINTS / 3 LINES\n\n"
    )

    print(
      summary(model_broken_3lines)
    )

    cat(
      "\n\nMODEL COMPARISON\n\n"
    )

    print(
      model_comparison,
      row.names = FALSE
    )

    cat(
      "\n\nTWO-LINE EQUATIONS\n\n"
    )

    print(
      parameter_table_2lines,
      row.names = FALSE
    )

    print(
      breakpoint_table_2lines,
      row.names = FALSE
    )

    cat(
      "\n\nTHREE-LINE EQUATIONS\n\n"
    )

    print(
      parameter_table_3lines,
      row.names = FALSE
    )

    print(
      breakpoint_table_3lines,
      row.names = FALSE
    )

    cat(
      sprintf(
        paste0(
          "\nAICc(3 lines) - AICc(2 lines) = %.6f\n"
        ),
        delta_AICc_3_vs_2
      )
    )
  },
  file = file.path(
    output_dir,
    "model_results_3_models.txt"
  )
)


cat(
  "\n\nAnalysis completed. Results saved in:\n",
  output_dir,
  "\n"
)
