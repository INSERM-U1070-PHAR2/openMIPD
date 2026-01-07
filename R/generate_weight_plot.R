#' Plot to show average model weights for categories/quantiles of a covariate in the test data
#'
#' @param data The generated test dataset with the columns MODEL, WEIGHT and covariates
#' @param cov Covariate
#' @param n_quantiles If continuous covariate, to how many quantiles should it be divided to display weights. Default is 4.
#' @param is_categorical Is the covariate categorical? TRUE or FALSE, default is FALSE
#' @returns A plot showing avareage weights for each model in the covariates categories/quantiles
#' @export
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @examples
#' weight_BURN <- generate_weight_plot(test_data, "BURN", is_categorical = TRUE)

generate_weight_plot <- function(data, cov, n_quantiles = 4, is_categorical = FALSE) {

  data_clean <- data |>
    filter(!is.na(.data[[cov]]) & !is.na(WEIGHT) & !is.na(MODEL))

  if (!is_categorical) {
    # Create quantiles for continuous coviables
    data_clean <- data_clean |>
      mutate(quantile_group = ntile(.data[[cov]], n_quantiles))

    # Calculate quantile ranges
    quantile_ranges <- quantile(data_clean[[cov]], probs = seq(0, 1, length.out = n_quantiles + 1), na.rm = TRUE)
    quantile_labels <- paste0(
      round(quantile_ranges[-(n_quantiles + 1)], 2), "-",
      round(quantile_ranges[-1], 3)
    )

    # Add quantile labels
    data_clean <- data_clean |>
      mutate(quantile_label = factor(quantile_group, levels = 1:n_quantiles, labels = quantile_labels))

    x_cov <- "quantile_label"
    x_label <- paste0(cov, " Quantile")
  } else {
    # Convert categorical coviables to factors
    data_clean <- data_clean |>
      mutate(!!cov := factor(.data[[cov]]))

    x_cov <- cov
    x_label <- cov
  }

  # Calculate average yeight
  average_weights <- data_clean |>
    group_by(MODEL, .data[[x_cov]]) |>
    summarise(avg_WEIGHT = mean(WEIGHT, na.rm = TRUE), .groups = "drop")

  # Bar plot
  plot <- ggplot(average_weights, aes(x = .data[[x_cov]], y = avg_WEIGHT, fill = MODEL)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(
      title = paste("Average model weight for", cov, ifelse(is_categorical, "categories", "quantiles")),
      x = x_label,
      y = "Average weight",
      fill = "MODEL"
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    ylim(0, 1) +
    theme(
      legend.position = "bottom",
      text = element_text(size = 12),
      axis.title = element_text(size = 18),
      plot.title = element_text(size=22)
    )

  return(plot)
}
