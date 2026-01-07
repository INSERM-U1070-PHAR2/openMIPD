#' Calculate relative root mean squared error and relative mean prediction error for the test data results.
#'
#' @param target_variable A secondary PK parameter (CMAX - maximal concentration, AUC - exposure or CMIN - trough concentration, ft_above_MIC - time above minimal inhibitory concentration) or CONCENTRATION.
#' @param test_data test data containing the columns ID, WEIGHTED_PREDICTION, and depending on the target variable CMAX/CMIN/AUC or TIME and DV for CONCENTRATION.
#' @returns A tibble with rRMSE and MPE.
#' @export
#' @import dplyr
#' @examples
#' compute_performance_metrics(test_data, target_variable = "CMAX")

compute_performance_metrics <- function(target_variable, test_data) {
  # Calculate mean of the target variable
  mean_target <- test_data |>
    summarise(mean_target = mean(.data[[target_variable]], na.rm = TRUE)) |>
    pull(mean_target)

  # Conditional grouping based on the target variable
  if (target_variable %in% c("CMAX", "AUC", "CMIN", "ft_above_MIC")) {
    target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

    # Group by ID and calculate performance metrics
    performance_metrics <- test_data |>
      group_by(ID) |>
      summarise(
        # Relative Root Mean Squared Error (RMSE)
        RMSE = sqrt(mean((.data[["WEIGHTED_PREDICTION"]] - .data[[target_column]])^2, na.rm = TRUE, finite = TRUE)) / mean(.data[[target_column]], na.rm = TRUE),

        # Relative Mean Prediction Error (MPE)
        MPE = mean(abs((.data[["WEIGHTED_PREDICTION"]] - .data[[target_column]]) / .data[[target_column]]), na.rm = TRUE, finite = TRUE)
      ) |>
      ungroup()
  } else if (target_variable == "CONCENTRATION") {
    # Group by ID and TIME
    performance_metrics <- test_data |>
      group_by(ID, TIME) |>
      summarise(
        # Relative Root Mean Squared Error (RMSE)
        RMSE = sqrt(mean((.data[["WEIGHTED_PREDICTION"]] - .data[["DV"]])^2, na.rm = TRUE, finite = TRUE)) / mean_target,

        # Relative Mean Prediction Error (MPE)
        MPE = mean(abs((.data[["WEIGHTED_PREDICTION"]] - .data[["DV"]]) / .data[["DV"]]), na.rm = TRUE, finite = TRUE)
      ) |>
      ungroup()
  } else {
    stop("Invalid target variable specified.")
  }

  # Remove infinite MPE values
  performance_metrics <- performance_metrics |>
    mutate(MPE = ifelse(is.infinite(MPE), NA, MPE))

  # Average the performance metrics
  final_metrics <- performance_metrics |>
    summarise(
      RMSE_avg = mean(RMSE, na.rm = TRUE, finite = TRUE),
      MPE_avg = mean(MPE, na.rm = TRUE, finite = TRUE)
    )

  return(final_metrics)
}

