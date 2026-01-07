#' Test classification tree model ensembling to attribute patient individualized model weights.
#'
#' @param test_data Test dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX), predicted target variable (for example CMAX_PRED), TIME if the target_variable is CONCENTRATION.
#' @param train_results Results of the algorithm training.
#'
#' @returns The test data with the attributed model weights.
#' @export
#' @import tidyr 
#' @import dplyr 
#' @import probably 
#' @import rlang
#' @examples
#' results <- classification_tree_model_ensembling_test(test_data = AMOX_CMAX_TEST, train_results = train_results)

classification_tree_model_ensembling_test <- function(test_data, train_results) {

  rpart_list <- train_results$cfit_list
  leaf_results <- train_results$leaf_results
  continuous_cov <- train_results$continuous_cov
  categorical_cov <- train_results$categorical_cov
  target_variable <- train_results$target_variable

  # Initialize score column
  test_data <- test_data |> mutate(score = NA) |> mutate(across(categorical_cov,~as.factor(.x)))

  target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  # Identify categorical covariates with more than two categories as they need to be treated like continuous covariates
  multi_category_covs <- categorical_cov[sapply(test_data[categorical_cov], function(x) length(unique(x)) > 2)]
  categorical_cov <- setdiff(categorical_cov, multi_category_covs)
  continuous_cov <- c(continuous_cov, multi_category_covs)
  models <- unique(test_data$MODEL)
  test_data_pooled <- NULL

  for(model in models){
    cfit <- rpart_list[[model]]
    test_data_filtered <-dplyr::filter(test_data,MODEL == model)
    test_data_filtered$score <- predict(cfit,test_data_filtered) |>
      as_tibble() |> pull(YES)
    test_data_pooled <- bind_rows(test_data_pooled,test_data_filtered)
  }

  test_data_pooled <- test_data_pooled |> arrange(ID,MODEL)
  # Determine grouping based on target variable
  group_vars <- if (target_variable == "CONCENTRATION") c("ID", "TIME") else "ID"

  # Standardize scores within each group
  test_data <- test_data_pooled |>
    group_by(across(all_of(group_vars))) |>
    mutate(
      score_sum = sum(score, na.rm = TRUE),
      WEIGHT = score / score_sum
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

