#' k nearest neighbors (KNN) for MIPD
#'
#' @param train A preprocessed (ml_data_preprocess function) training dataset with one line per patient containing covariates, the administered dose, and interdose interval.
#' @param continuous_cov List of continuous covariates.
#' @param categorical_cov List of categorical covariates.
#' @param seed Random seed.
#' @param grid_size Grid size for hyperparameter tuning.
#' @returns Trained KNN model, tuning autoplots, and tuned hyperparameters.
#' @export
#' @import tidymodels
#' @import vip
#' @import shapviz
#' @import doParallel
#' @import ggplot2
#' @import brulee
#' @import tidyr
#' @import purrr
#' @import tidyverse
#' @import parallel
#' @import kknn
#' @import stacks
#' @import yardstick
#' @import recipes
#' @import parsnip
#' @import rsample
#' @examples
#' results <- knn_train(train = AMOX_CMIN_TRAIN, continuous_cov = c("WT", "CRCL"), categorical_cov = c("BURN", "OBESE"))

knn_train <- function(train,
                      continuous_cov,
                      categorical_cov,
                      grid_size = 10,
                      seed = 1991) {

  set.seed(seed)

  reg_metrics <- metric_set(mae)

  # predictors (covariates + number of daily administrations)
  predictors <- c(continuous_cov, categorical_cov, "INF")
  formula <- as.formula(paste("log_DOSE_TARGET ~", paste(predictors, collapse = " + ")))

  data_recipe <- recipe(formula, data = train) |>
    step_mutate_at(all_of(categorical_cov), fn = ~factor(.)) |>
    step_dummy(all_nominal_predictors()) |>
    step_normalize(all_numeric_predictors())

  # Model specifications
  knn_spec <- nearest_neighbor(mode = "regression", neighbors = tune()) |>
    set_engine("kknn")

  # Cross-validation
  set.seed(seed)
  folds <- vfold_cv(train, v = 10, strata = strata_dosing_regimen)

  # Parallel
  num_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(num_cores, type = "PSOCK",
                              rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
  doParallel::registerDoParallel(cl)

  # Helper that tunes a model and returns useful objects
  tune_predict <- function(spec, name, recipe_obj, grid_size) {
    wf <- workflow() %>%
      add_model(spec) %>%
      add_recipe(recipe_obj)

    tune_res <- tune_grid(
      wf,
      resamples = folds,
      grid = grid_size,
      metrics = reg_metrics,
      control = control_stack_grid(),
    )

    best_params <- select_best(tune_res, metric = "mae") # tuning metric: mean average error
    final_wf <- finalize_workflow(wf, best_params)
    final_fit <- fit(final_wf, data = train)

    list(
      name = name,
      workflow = wf,
      recipe = recipe_obj,
      tune_res = tune_res,
      best_params = best_params,
      final_wf = final_wf,
      final_fit = final_fit
    )
  }

  # Tune all models
  tune_knn <- tune_predict(knn_spec, "KNN", data_recipe, grid_size)

  # Stop parallel processing
  stopCluster(cl)
  registerDoSEQ()

  # Tuning autoplots
  tuning_knn <- autoplot(tune_knn$tune_res, metric = "mae", scientific = FALSE) +
    theme_bw() +
    ggtitle("Hyperparameter tuning - KNN")

  # Final workflows
  final_knn_wf <- tune_knn$final_wf

  # Final fitted workflows
  final_knn_fit <- tune_knn$final_fit

  # Results
  return(list(tune_plot_knn = tuning_knn, final_wf_knn = final_knn_wf, tune_res_knn = tune_knn$tune_res, final_knn_fit = final_knn_fit))

}
