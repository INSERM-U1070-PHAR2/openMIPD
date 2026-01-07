#' Factor Analysis of Mixed Data (FAMD) test to attribute model weights based on the similiarty of a patient's covariates to the model development cohorts
#'
#' @param test A test dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX), predicted target variable (for example CMAX_PRED), MODEL name used for the simulation.
#' @param train_results The FAMD model resulting from the famd_train function.
#' @returns A goodness of fit (GOF) plot and the test dataset with the attributed weights and ensembled concentrations.
#' @export
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import FNN
#' @import factoextra
#' @import here
#' @import probably
#' @import RSpectra
#' @import purrr
#' @import stringr
#' @import FactoMineR
#'
#' @examples
#' result<- famd_train(test = test, train_results = train_results)

famd_test <- function(test, train_results) {

  categorical_cov <- train_results$categorical_cov
  continuous_cov <- train_results$continuous_cov
  target_variable <- train_results$target_variable
  pca_result <- train_results$pca_result
  optimal_components <- train_results$optimal_components
  centroids <- train_results$centroids
  cov_matrices <- train_results$cov_matrices

  # Define covariates
  covariates <- c(continuous_cov, categorical_cov)

  # Get the different models
  model_groups <- unique(test$MODEL)

  target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  test <- test %>%
    mutate(across(all_of(categorical_cov), as.factor))

  test_latent <- predict(pca_result, newdata = test[, covariates])$coord[, 1:optimal_components]

  # Compute Mahalanobis distances and weights
  weights_list <- apply(test_latent, 1, function(new_patient_pca) {
    distances <- sapply(names(centroids), function(model) {
      sqrt(mahalanobis(new_patient_pca, centroids[[model]], cov_matrices[[model]]))
    })
    weights <- 1 / distances
    weights / sum(weights)
  })

  # Convert weights to a data frame
  weights_df <- as.data.frame(t(weights_list))
  weights_df <- distinct(weights_df)

  # Bind weights with relevant rows of test_data
  test_data <- test |>
    dplyr::select(ID, all_of(covariates)) |>
    distinct()
  test_weights <- cbind(test_data, weights_df)

  # Change to long format to create model WEIGHT column
  test_weights_long <- test_weights |>
    pivot_longer(cols = all_of(model_groups), names_to = "MODEL", values_to = "WEIGHT")

  # Merge with the original test data
  test_data <- merge(test, test_weights_long,
                   by = c("ID", all_of(covariates), "MODEL"),
                   all.x = TRUE)

  # Get final ensembled prediction
  prediction_col <- paste0(target_variable, "_PRED")
  weighted_prediction_col <- paste0("WEIGHTED_PREDICTION")

  # Determine how to group data (use TIME for concentration)
  group_vars <- if (target_variable == "CONCENTRATION") c("ID", "TIME") else "ID"

  test_data <- test_data |>
    mutate(WEIGHTED_VALUE = .data[[prediction_col]] * WEIGHT) |>
    group_by(!!!syms(group_vars)) |>
    mutate("{weighted_prediction_col}" := sum(WEIGHTED_VALUE)) |>
    ungroup()

  # Goodness of fit plot
  GOF_plot <- cal_plot_regression(test_data, !!sym(weighted_prediction_col), truth = !!sym(target_column))

  return(list(test_results = test_data, GOF_plot = GOF_plot))

}
