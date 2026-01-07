#' Test regression tree model ensembling to attribute subjects individualized model weights
#'
#' @param train_results Results of the classification tree training.
#' @param test_data Test dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX_IND), predicted target variable (for example CMAX_PRED), TIME if the target_variable is CONCENTRATION.
#' @returns The test data with the attributed model weights.
#' @export
#' @import tidyr
#' @import dplyr
#' @import probably
#' @import rlang
#' @examples
#' results <- regression_tree_model_ensembling_test(test_data = AMOX_CMAX_TEST, train_results = train_results)


regression_tree_model_ensembling_test <- function(test_data, train_results) {

    model_trees <- train_results$model_trees
  categorical_cov <- train_results$categorical_cov
  continuous_cov <- train_results$continuous_cov
  target_variable <- train_results$target_variable

  target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  models <- unique(test_data$MODEL)

  # Predict RATIO for test_data using the correct model tree
  test_data$RATIO <- NA  # Initialize column

  for (model in models) {
    if (!model %in% names(model_trees)) {
      warning(paste("No trained tree found for model:", model))
      next
    }
    # Filter test_data for the current model
    test_subset <- test_data %>% dplyr::filter(MODEL == model)

      # Get the corresponding trained model
    cfit <- model_trees[[model]]

      # Convert categorical variables to factors
    for (cat_var in categorical_cov) {
      test_subset[[cat_var]] <- as.factor(test_subset[[cat_var]])
    }

      # Predict RATIO and add to original test data
    test_data$RATIO[test_data$MODEL == model] <- predict(cfit, newdata = test_subset)
  }

  # Determine grouping based on target variable
  group_vars <- if (target_variable == "CONCENTRATION") c("ID", "TIME") else "ID"

  test_data <- test_data |>
    mutate(RATIO = ifelse(RATIO == 0, 0.00001, RATIO)) |>  # Replace 0 with a small value to be able to take the reciprocal
    mutate(RATIO_REC = abs(1 / RATIO))

  # Standardize scores within each group
  test_data <- test_data |>
    group_by(across(all_of(group_vars))) |>
    mutate(
      WEIGHT = RATIO_REC / sum(RATIO_REC, na.rm = TRUE)
    ) |>
    ungroup()

  # Calculate the weighted prediction
  pred_col <- paste0("WEIGHTED_", target_variable)

  test_data <- test_data |>
    mutate(!!pred_col := WEIGHT * .data[[paste0(target_variable, "_PRED")]])

  if (!pred_col %in% names(test_data)) {
    stop(paste("Column", pred_col, "not found in test_data after mutation."))
  }

  # Compute final weighted prediction
  test_data <- test_data |>
    group_by(across(all_of(group_vars))) |>
    mutate(WEIGHTED_PREDICTION = sum(.data[[pred_col]], na.rm = TRUE)) |>
    ungroup()

  # Goodness of fit plot
  GOF_plot <- cal_plot_regression(test_data, !!sym("WEIGHTED_PREDICTION"), truth = !!sym(target_column))

  # Boxplot for predicted vs observed ratio with the bioequivalence range
  Boxplot <- ggplot(test_data, aes(y = !!sym("WEIGHTED_PREDICTION") / !!sym(target_column) * 100)) +
    geom_boxplot() +
    scale_y_log10() +
    geom_hline(yintercept = c(80, 125), color = "red", linetype = "dashed") +
    labs(y = "ratio in %", title = "Predicted/observed ratio") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1, size = 16),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size=24),
    )

  return(list(test_results = test_data, GOF_plot = GOF_plot, Boxplot = Boxplot))
}
