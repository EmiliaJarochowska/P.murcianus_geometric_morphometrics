# Function to calculate group-wise slopes
calculate_slopes <- function(data, x_var, y_var, group_var) {
  groups <- unique(data[[group_var]])
  cat("Slopes by", group_var, ":\n")
  
  for (grp in groups) {
    if (!is.na(grp)) {
      sub_data <- dplyr::filter(data, !!rlang::sym(group_var) == grp)
      if (nrow(sub_data) > 1) {
        model <- stats::lm(formula(paste(y_var, "~", x_var)), data = sub_data)
        slope <- stats::coef(model)[x_var]
        r2 <- summary(model)$r.squared
        cat("  ", group_var, ":", grp, 
            "| Slope =", round(slope, 4), 
            "| R² =", round(r2, 3), "\n")
      }
    }
  }
}