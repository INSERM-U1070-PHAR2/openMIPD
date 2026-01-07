#' Train weighed model ensembling to calculate model and subgroup influence scores
#'
#' @param target_variable A secondary PK parameter (CMAX - maximal concentration, AUC - exposure or CMIN - trough concentration), time above minimal inhibitory concentration (ft_above_MIC) or CONCENTRATION.
#' @param data A dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX), predicted target variable (for example CMAX_PRED), and MODEL name used for the simulation.
#' @param pen_RMSE - Negative penalty for relative root mean squared error (RMSE), default = -25.
#' @param pen_MPE - Negative penalty for relative mean prediction error (MPE), default = -3.
#' @param subgroup_inf_lower_bound Lower bound of the percentage interval within which a prediction is considered correct for subgroup influence calculation (Defaults to 80 %).
#' @param subgroup_inf_upper_bound Upper bound of the percentage interval within which a prediction is considered correct for subgroup influence calculation (Defaults to 125 %).
#' @param continuous_cov = List of continuous covariates.
#' @param categorical_cov = List of categorical covariates.
#' @returns A table with RMSE, MPE, model influence and subgroup influence for each model and covariate quantile/category.
#' @export
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @examples
#' weighed_model_ensembling_train(data = data, target_variable = "CMAX", pen_RMSE = -25, pen_MPE = -3, continuous_cov = c("WT", "CRCL"), categorical_cov = c("SEX", "BURN"))


weighed_model_ensembling_train <- function(data, target_variable, pen_RMSE = -25, pen_MPE = -3,
                                           subgroup_inf_lower_bound = 80,
                                           subgroup_inf_upper_bound = 125,
                                           continuous_cov, categorical_cov) {

  target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  # Function to assign quantiles or categories
  assign_quantiles_categories <- function(data, covariate, is_categorical, levels = NULL) {
    if (is_categorical) {
      data[[paste0("QUANTILE_", covariate)]] <- as.character(data[[covariate]])
    } else {
      quantile_bins <- quantile(data[[covariate]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
      data[[paste0("QUANTILE_", covariate)]] <- cut(
        data[[covariate]], breaks = quantile_bins, include.lowest = TRUE,
        labels = paste0(signif(head(quantile_bins, -1), 3), "-", signif(tail(quantile_bins, -1), 3))
      )
    }
    return(data)
  }

  # Assign quantiles/categories
  for (var in continuous_cov) {
    data <- assign_quantiles_categories(data, var, FALSE)
  }
  for (var in categorical_cov) {
    data <- assign_quantiles_categories(data, var, TRUE)
  }

  # Calculate differences and ratios
  if (target_variable %in% c("CMAX", "CMIN", "AUC", "ft_above_MIC")) {
    differences <- data |>
      group_by(ID) |>
      mutate(diff = .data[[paste0(target_variable, "_PRED")]] - .data[[paste0(target_column)]]) |>
      mutate(ratio = (.data[[paste0(target_variable, "_PRED")]] /.data[[paste0(target_column)]])*100) |>
      ungroup()
  } else if (target_variable == "CONCENTRATION") {
    differences <- data |>
      group_by(ID, TIME) |>
      mutate(diff = .data[[paste0(target_variable, "_PRED")]] - .data[[paste0(target_column)]]) |>
      mutate(ratio = (.data[[paste0(target_variable, "_PRED")]] /.data[[paste0(target_column)]])*100) |>
      ungroup()
  }

  # Calculate subgroup influence
  calculate_proportions <- function(differences, quantile_col,
                                    subgroup_inf_lower_bound = 80,
                                    subgroup_inf_upper_bound = 125) {
    differences |>
      dplyr::filter(!is.infinite(ratio) & !is.na(ratio)) |>
      group_by(.data[[quantile_col]]) |>
      summarise(SUBGROUP_INFLUENCE = sum(ratio < subgroup_inf_lower_bound | ratio > subgroup_inf_upper_bound, na.rm = TRUE) / n(), .groups = "drop") |>
      dplyr::mutate(COVARIATE = gsub("QUANTILE_", "", quantile_col))
  }

  proportion_data <- bind_rows(
    lapply(c(paste0("QUANTILE_", continuous_cov), paste0("QUANTILE_", categorical_cov)),
           function(var) calculate_proportions(differences, var))
  )

  proportion_data <- proportion_data |>
    pivot_longer(cols = starts_with("QUANTILE"), names_to = "QUANTILE", values_to = "Quantile") |>
    dplyr::filter(!is.na(Quantile)) |>
    dplyr::select(COVARIATE, QUANTILE = Quantile, SUBGROUP_INFLUENCE)

  # Calculate model scores
  calculate_scores <- function(data, quantile_col) {
    scores <- data |>
      group_by(MODEL, .data[[quantile_col]]) |>
      summarise(
        MPE = mean(abs(diff) / get(target_column), na.rm = TRUE),
        RMSE = sqrt(mean(diff^2, na.rm = TRUE)),
        mean_value = mean(get(target_column), na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        RMSE = RMSE / mean_value,
        score_MPE = exp(abs(MPE) * pen_MPE),
        score_RMSE = exp(abs(RMSE) * pen_RMSE),
        MODEL_INFLUENCE = score_RMSE * score_MPE,
        COVARIATE = gsub("QUANTILE_", "", quantile_col)
      )
    return(scores)
  }

  model_data <- bind_rows(
    lapply(c(paste0("QUANTILE_", continuous_cov), paste0("QUANTILE_", categorical_cov)),
           function(var) calculate_scores(differences, var))
  )

  model_data <- model_data |>
    pivot_longer(cols = starts_with("QUANTILE"), names_to = "QUANTILE", values_to = "Quantile") |>
    dplyr::filter(!is.na(Quantile)) |>
    dplyr::select(MODEL, MPE, RMSE, COVARIATE, QUANTILE = Quantile, MODEL_INFLUENCE)

  # Merge proportion and model data
  all_scores <- merge(proportion_data, model_data, by = c("COVARIATE", "QUANTILE"), all.x = TRUE) |>
    distinct()

  # For continuous covariates, change upper limit to infinity and lower limit to 0 to be able to handle new data which is out of the training range
  all_scores <- all_scores %>%
    mutate(
      LABEL = paste0(COVARIATE, "_", as.integer(factor(QUANTILE)))
    ) %>%
    mutate(
      lower = if_else(
        COVARIATE %in% continuous_cov,
        as.numeric(sub("-.*", "", QUANTILE)),
        NA_real_
      ),
      upper = if_else(
        COVARIATE %in% continuous_cov,
        as.numeric(sub(".*-", "", QUANTILE)),
        NA_real_
      )
    ) %>%
    group_by(COVARIATE) %>%
    mutate(
      upper = if_else(
        COVARIATE %in% continuous_cov & row_number() == n(),
        Inf, upper
      ),
      lower = if_else(
        COVARIATE %in% continuous_cov & row_number() == 1,
        0, lower
      )
    ) %>%
    ungroup()

  # Create a summary plot
  MPE_plot <- ggplot(all_scores, aes(x = MODEL, y = MPE, fill = MODEL)) +
    geom_col(position = "dodge") +
    facet_wrap(~ LABEL, scales = "free_y") +
    labs(title = paste("MPE for each model & subgroup (", target_variable, ")")) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), strip.text = element_text(size = 16))

  RMSE_plot <- ggplot(all_scores, aes(x = MODEL, y = RMSE, fill = MODEL)) +
    geom_col(position = "dodge") +
    facet_wrap(~ LABEL, scales = "free_y") +
    labs(title = paste("RMSE for each model & subgroup (", target_variable, ")")) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), strip.text = element_text(size = 16))

  # Plot for subgroup influences
  all_scores_subgroup <- all_scores %>%
    dplyr::select(SUBGROUP_INFLUENCE, LABEL) %>%
    distinct()

  SUBGROUP_INFLUENCE_plot <- ggplot(all_scores_subgroup, aes(x = LABEL, y = SUBGROUP_INFLUENCE)) +
    geom_col(fill = "steelblue", position = "dodge") +
    labs(
      x = "Subgroup",
      y = "Subgroup influence",
      title = "Subgroup influences"
    ) +
    ylim(0,1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 24),
      axis.text = element_text(size = 12)
    )

  return(list(all_scores = all_scores, MPE_plot = MPE_plot, RMSE_plot = RMSE_plot, SUBGROUP_INFLUENCE_plot = SUBGROUP_INFLUENCE_plot, continuous_cov = continuous_cov, categorical_cov = categorical_cov, target_variable = target_variable))
}

