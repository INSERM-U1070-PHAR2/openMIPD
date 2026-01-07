#' Support Vector Machine (SVM) for MIPD
#'
#' @param train A preprocessed (ml_data_preprocess function) training dataset with one line per patient containing covariates, the administered dose, and interdose interval.
#' @param continuous_cov List of continuous covariates.
#' @param categorical_cov List of categorical covariates.
#' @param seed Random seed.
#' @param grid_size Grid size for hyperparameter tuning.
#' @returns Tuning autoplots, tuned hyperparameters, and variable importance plot, and trained model.
#' @export
#' @import yardstick
#' @import recipes
#' @import parsnip
#' @import rsample
#' @import workflows

#' @examples
#' results <- svm(train = AMOX_CMIN_TRAIN, continuous_cov = c("WT", "CRCL"), categorical_cov = c("BURN", "OBESE"))

svm_train <- function(train,
                continuous_cov,
                categorical_cov,
                grid_size = 10,
                seed = 1991) {

  set.seed(seed)

  library(tidymodels)
  library(vip)
  library(shapviz)
  library(doParallel)
  library(ggplot2)
  library(brulee)
  library(tidyr)
  library(purrr)
  library(tidyverse)
  library(parallel)
  library(stacks)

  reg_metrics <- metric_set(mae)

  # predictors (covariates + number of daily administrations)
  predictors <- c(continuous_cov, categorical_cov, "INF")
  formula <- as.formula(paste("log_DOSE_TARGET ~", paste(predictors, collapse = " + ")))

  data_recipe <- recipe(formula, data = train) |>
    step_mutate_at(all_of(categorical_cov), fn = ~factor(.)) |>
    step_dummy(all_nominal_predictors()) |>
    step_normalize(all_numeric_predictors())

  # Model specifications
  svm_spec <- svm_rbf(mode = "regression", cost = tune(), rbf_sigma = tune()) |>
    set_engine("kernlab")

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

    best_params <- select_best(tune_res, metric = "mae")
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
  tune_svm <- tune_predict(svm_spec, "SVM", data_recipe, grid_size)

  # Stop parallel processing
  stopCluster(cl)
  registerDoSEQ()

  # Tuning autoplots
  tuning_svm <- autoplot(tune_svm$tune_res, metric = "mae", scientific = FALSE) +
    theme_bw() +
    ggtitle("Hyperparameter tuning - SVM")

  # Final workflows
  final_svm_wf <- tune_svm$final_wf

  # Final fitted workflows
  final_svm_fit <- tune_svm$final_fit

  # Results
  return(list(tune_plot_svm = tuning_svm, final_wf_svm = final_svm_wf, tune_res_svm = tune_svm$tune_res, final_svm_fit = final_svm_fit))

}
