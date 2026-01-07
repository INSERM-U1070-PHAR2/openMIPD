#' Test weighed model ensembling to attribute subjects individualized model weights
#'
#' @param train_results Results of the weighed model ensembling training.
#' @param test_data Test dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX), predicted target variable (for example CMAX_PRED), TIME if the target_variable is CONCENTRATION.
#' @returns The test data with the attributed model weights and ensembled concentrations, goodnes of fit (GOF) plot and predicted/observed concentration ratio plot.
#' @export
#' @import probably
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import rlang
#' @examples
#' result <- weighed_model_ensembling_test(train_results = train_results, test_data = AMOX_CMAX_TEST)


# Function to compute model influences and weighted predictions
weighed_model_ensembling_test <- function(test_data, train_results) {

  all_scores <- train_results$all_scores
  categorical_cov <- train_results$categorical_cov
  continuous_cov <- train_results$continuous_cov
  target_variable <- train_results$target_variable

  target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  get_subgroup_influence <- function(covariate_value, covariate_name, all_scores) {
    covariate_scores <- all_scores |> dplyr::filter(COVARIATE == covariate_name)
    match <- covariate_scores |> dplyr::filter(covariate_value >= lower & covariate_value <= upper)

    if (nrow(match) > 0) {
      return(first(match$SUBGROUP_INFLUENCE))
    } else {
      return(NA)
    }
  }

  get_model_influence_for_range <- function(value, covariate, model, all_scores) {
    model_influence <- all_scores |>
      dplyr::filter(MODEL == model & COVARIATE == covariate & value >= lower & value <= upper) |>
      dplyr::select(MODEL_INFLUENCE) |>
      pull()
    if (length(model_influence) == 1) return(model_influence) else return(NA)
  }

  # Put together all covariates
  all_covariates <- c(continuous_cov, categorical_cov)


  get_subgroup_influence_new <- function(data, score_data, covariates) {

    convert_quantiles_to_breaks <- function(quantiles) {
      if(sum(grepl("-",quantiles)) >=1){
        str_split(quantiles, pattern = "-") |> unlist() |> as.numeric()
      } else {
        c(-Inf,quantiles |> as.numeric())
      }
    }

    score_data_reduced <- all_scores |> select(COVARIATE,QUANTILE,SUBGROUP_INFLUENCE) |>
      distinct() |>
      pivot_wider(names_from = COVARIATE,
                  names_glue = "{COVARIATE}_QUANTILE",
                  values_from = QUANTILE)

    data |>
      mutate(across(
        all_of(all_covariates),
        .fns = ~ cut(.x, breaks = convert_quantiles_to_breaks((
          all_scores$QUANTILE[all_scores$COVARIATE == cur_column()]
        )) |> unique() |> sort(),
        labels = unique(all_scores$QUANTILE[all_scores$COVARIATE == cur_column()]))
        , .names = "{.col}_QUANTILE"
      )) |>
      mutate(across(
        all_of(paste0(all_covariates,"_QUANTILE")),
        .fns = ~ score_data_reduced$SUBGROUP_INFLUENCE[match(.x,score_data_reduced[[cur_column()]])]
        , .names = "{.col}_SUBGROUP_INFLUENCE"
      )) |>
      rename_with(.fn = ~gsub("_QUANTILE_SUBGROUP_INFLUENCE","_SUBGROUP_INFLUENCE",.x) )
  }

  get_model_influence_for_range_new <- function(data,score_data,covariates) {
    score_data_reduced <- all_scores |> select(COVARIATE,QUANTILE,MODEL,MODEL_INFLUENCE) |>
      distinct() |>
      mutate(MODEL_QUANTILE = paste0(MODEL,"-",QUANTILE)) |>
      pivot_wider(names_from = COVARIATE,
                  names_glue = "{COVARIATE}_MODEL_QUANTILE",
                  values_from = MODEL_QUANTILE)

    data |>
      mutate(across(
        all_of(paste0(all_covariates,"_QUANTILE")),
        .fns = ~ paste0(.data$MODEL,"-",.x)
        , .names = "{.col}_MODEL_QUANTILE")) |>
      rename_with(.fn = ~gsub("_QUANTILE_MODEL_QUANTILE","_MODEL_QUANTILE",.x) ) |>
      mutate(across(
        all_of(paste0(all_covariates,"_MODEL_QUANTILE")),
        .fns = ~ score_data_reduced$MODEL_INFLUENCE[match(.x,score_data_reduced[[cur_column()]])]
        , .names = "{.col}_MODEL_INFLUENCE"
      )) |>
      rename_with(.fn = ~gsub("_MODEL_QUANTILE_MODEL_INFLUENCE","_MODEL_INFLUENCE",.x) )
  }

  test_data <- get_subgroup_influence_new(data = test_data,
                                          score_data = all_scores,
                                          covariates = all_covariates) |>
    get_model_influence_for_range_new(score_data = all_scores,
                                      covariates = all_covariates)

  compute_colname <- function(colname,suffix){paste0(colname,suffix)}

  weight_table <- test_data |>
    mutate(across(all_of(all_covariates),
                  .fns = ~ get(paste0(cur_column(),"_MODEL_INFLUENCE")) * get(paste0(cur_column(),"_SUBGROUP_INFLUENCE")),
                  .names = "{.col}_MODEL_INFLUENCE_FINAL")
    )  |>
    group_by(ID,MODEL) |>
    select((ends_with("_MODEL_INFLUENCE_FINAL"))) |>
    pivot_longer(cols = ends_with("_MODEL_INFLUENCE_FINAL")) |>
    mutate(
      WEIGHT = sum(value)) |>
    ungroup() |>
    select(ID,MODEL,WEIGHT) |>
    distinct()

  # Compute weight as the sum of all model influences
  test_data2 <- test_data |>
    left_join(weight_table)

  # Use also time for grouping if the target variable is concentration
  if (target_variable == "CONCENTRATION") {
    test_data2 <- test_data2 |> group_by(ID, TIME)
  } else {
    test_data2 <- test_data2 |> group_by(ID)
  }

  test_data3 <- test_data2 |>
    mutate(
      WEIGHT = WEIGHT / sum(WEIGHT, na.rm = TRUE)
    ) |>
    ungroup()

  # Calculate weighted prediction
  prediction_col <- paste0(target_variable, "_PRED")  # Dynamically get the column name
  weighted_prediction_col <- paste0("WEIGHTED_PREDICTION")

  # Ensemble the predictions
  group_vars <- if (target_variable == "CONCENTRATION") c("ID", "TIME") else "ID"
  test_data <- test_data3 |>
    mutate(WEIGHTED_VALUE = !!sym(prediction_col) * WEIGHT) |>
    group_by(across(all_of(group_vars))) |>
    mutate("{weighted_prediction_col}" := sum(WEIGHTED_VALUE)) |>
    ungroup()

  # Goodness of fit plot (predictions as a function of observations)
  GOF_plot <- cal_plot_regression(test_data, !!sym(weighted_prediction_col), truth = !!sym(target_column))

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

