#' Train regression tree informed model ensembling to calculate weights with log-transformed prediction/observation ratio as target variable
#'
#' @param target_variable A secondary PK parameter (CMAX - maximal concentration, AUC - exposure or CMIN - trough concentration, CONCENTRATION - complete concentration profileCONCEN).
#' @param data A dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX), predicted target variable (for example CMAX_PRED), MODEL name used for the simulation.
#' @param continuous_cov List of continuous covariates.
#' @param categorical_cov List of categorical covariates.
#' @param seed Set seed for reproductibility.
#' @param conservative_pruning If FALSE (default) take the highest CP with the lowest cross-validation error else take the highest CP whose cross-validation error is within 1 SD of the lowest cross-validation error.
#' @returns Trained regression trees.
#' @examples leaf_results <- regression_tree_model_ensembling_train(data = data, target_variable = "CMAX", continuous_cov = c("WT", "CRCL"), categorical_cov = c("SEX", "BURN"))
#' @import tidymodels
#' @import dplyr
#' @import rpart.plot
#' @export


regression_tree_model_ensembling_train <- function(data, target_variable, continuous_cov, categorical_cov, seed = 1991, conservative_pruning = FALSE) {

  set.seed(seed)

  target_column <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  # Function to calculate differences and ratios between predicted and observed values
  calculate_ratio <- function(data, target_column, target_variable) {
    ref_value <- data[[target_column]][1]
    data$RATIO <- log((data[[paste0(target_variable, "_PRED")]]) / ref_value)
    return(data)
  }

  if (target_variable %in% c("CMAX", "CMIN", "AUC", "ft_above_MIC")) {
    data <- data |>
      group_by(ID) |>
      do(calculate_ratio(., target_column, target_variable)) |>
      ungroup()
  } else if (target_variable == "CONCENTRATION") {
    data <- data |>
      group_by(ID, TIME) |>
      do(calculate_ratio(., "DV")) |>
      ungroup()
  }

  for (cat_var in categorical_cov) {
    data[[cat_var]] <- as.factor(data[[cat_var]])
  }

  models <- unique(data$MODEL)

  covariates <- c(continuous_cov, categorical_cov)

  leaf_results <- list()
  model_trees <- list()

  for (model in models) {
    sub <- data |>
     dplyr::filter(MODEL == model) |>
      select(all_of(covariates), RATIO)

    rec <- recipe(RATIO ~ ., data = sub)

    tree_spec <- decision_tree(
      cost_complexity = tune(),   # tune complexity parameter (=cp)
      tree_depth      = tune(),   # tune train depth
      min_n           = tune()    # tune minsplit/minbucket
    ) |>
      set_engine("rpart") |>
      set_mode("regression")

    wf <- workflow() |>
      add_model(tree_spec) |>
      add_recipe(rec)

    folds <- vfold_cv(sub, v = 10, strata = NULL)

    grid <- expand.grid(
      cost_complexity = seq(0.001, 0.02, by = 0.001),
      tree_depth      = c(1, 2, 3, 4, 5),
      min_n           = c(4)
    )

    res <- tune_grid(
      wf,
      resamples = folds,
      grid = grid,
      metrics = metric_set(rmse)
    )

    # Either take the highest CP with the lowest cross-validation error
    # Or be more conservative and take the highest CP whose cross-validation
    # error is within 1 SD of the lowest cross-validation error
    if (conservative_pruning == FALSE) {
      best_params <- select_best(res, metric = "rmse")
    } else {
      best_params <- select_by_one_std_err(res, metric = "rmse", cost_complexity, tree_depth, min_n)
    }

    final_wf <- finalize_workflow(wf, best_params)

    final_fit <- final_wf %>% fit(data = sub)

    # Extract rpart model for plotting
    fitted_rpart <- extract_fit_engine(final_fit)

    plot <- rpart.plot(fitted_rpart,
                       main = paste("Decision Tree for", model),
                       box.palette = "Blues"
    )

    model_trees[[model]] <- fitted_rpart
  }

  return(list(model_trees = model_trees, continuous_cov = continuous_cov, categorical_cov = categorical_cov, target_variable = target_variable))
}

