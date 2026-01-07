#' Four machine learning models are fit to make predictions for the test data - Support Vector Machine (SVM), k nearest neighbors (KNN), XGBoost, Random forest (RF)
#'
#' @param target_variable A secondary PK parameter (CMAX - maximal concentration, AUC - exposure or CMIN - trough concentration).
#' @param train A training dataset with one line per patient. The dataset should include the target variable, covariables and the DAILY dose and a column called II for interdose interval (set to 24 for continuous administration)
#' @param train A test dataset with one line per patient. The dataset should include the target variable, covariables and the DAILY dose.
#' @param continuous_cov List of continuous covariates.
#' @param categorical_cov List of categorical covariates.
#' @param target_concentration Target to reach (Cmax, Cmin or AUC).
#' @returns A test data with four additional columns containing the predictions for the four ML models.
#' @export
#' @import tidymodels
#' @import vip
#' @import shapviz
#' @import doParallel
#' @import ggplot2
#' @import ranger
#' @import kknn
#' @import brulee
#' @import stacks
#' @import tidyr
#' @import rpart.plot
#' @import purrr
#' @import stringr
#' @import tidyverse
#' @import yardstick
#' @import recipes
#' @import parsnip
#' @import rsample
#' @import workflows
#' @import tune
#' @import foreach
#' @examples
#' results <- machine_learning(train = AMOX_CMIN_TRAIN, test = AMOX_CMIN_TEST, continuous_cov = c("WT", "CRCL"), categorical_cov = c("BURN", "OBESE"), target_variable = "CMIN", target_concentration = 60)

machine_learning <- function(train,
                               test,
                               continuous_cov,
                               categorical_cov,
                               target_variable,
                               target_concentration) {

  set.seed(1991)

  target_variable <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

  reg_metrics <- metric_set(mae)

  # Preprocessing, convert administered dose to daily dose
  preprocess_data <- function(data) {
    data |>
      dplyr::mutate(
        INF = 24 / II,
        DOSE_ADM = DOSE_ADM * (24 / II),
        DOSE_TARGET = (target_concentration / .data[[target_variable]]) * DOSE_ADM,
        strata_dosing_regimen = as.factor(paste0(DOSE_ADM/INF,"mg q",FREQ,"h dur",DUR,"h")),
        log_DOSE_TARGET = log(DOSE_TARGET)
      )
  }

  train <- preprocess_data(train)
  test <- preprocess_data(test)

  # Recipe
  predictors <- c(continuous_cov, categorical_cov, "INF")
  formula <- as.formula(paste("log_DOSE_TARGET ~", paste(predictors, collapse = " + ")))

  # Recipe for SVM and KNN with data normalization
  data_recipe_svm_knn <- recipe(formula, data = train) |>
    step_mutate_at(all_of(categorical_cov), fn = ~factor(.)) |>
    step_dummy(all_nominal_predictors()) |>
    step_normalize(all_numeric_predictors())

  # Recipe for RF and XGB without normalizing the data
  data_recipe_rf_xgb <- recipe(formula, data = train) |>
    step_mutate_at(all_of(categorical_cov), fn = ~factor(.)) |>
    step_dummy(all_nominal_predictors())

  # Model specifications
  # Support Vector Machine
  svm_spec <- svm_rbf(mode = "regression", cost = tune(), rbf_sigma = tune()) |>
    set_engine("kernlab")

  # Random forest
  rf_spec <- rand_forest(mode = "regression", trees = tune(), mtry = tune()) |>
    set_engine("ranger", importance = "impurity")

  # XGBoost
  xgb_spec <- boost_tree(mode = "regression",
                         trees = tune(),
                         min_n = tune(),
                         tree_depth = tune(),
                         learn_rate = tune()) |>
    set_engine("xgboost")

  # k nearest neighbors
  knn_spec <- nearest_neighbor(mode = "regression", neighbors = tune()) |>
    set_engine("kknn")

  # Cross-validation (10-fold) stratified by dosing regimen
  set.seed(1995)
  folds <- vfold_cv(train, v = 10, strata = strata_dosing_regimen)

  # Parallel processing
  num_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(num_cores, outfile = "log.txt", type = "PSOCK",
                              rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
  doParallel::registerDoParallel(cl)

  # Make predictions after tuning
  tune_predict <- function(spec, name, resamples, grid_size = 10) {
    # select recipe based on model name
    chosen_recipe <- if (toupper(name) %in% c("RF", "XGB")) {
      data_recipe_rf_xgb
    } else {
      data_recipe_svm_knn
    }

    wf <- workflow() %>%
      add_model(spec) %>%
      add_recipe(chosen_recipe)

    tuned <- tune_grid(
      wf,
      resamples = folds,
      grid = grid_size,
      metrics = reg_metrics,   # <- add this line
      control = control_grid(save_pred = TRUE, verbose = FALSE)
    )

    # pick best params on mean absolute error
    best_params <- select_best(tuned, metric = "mae")

    final_wf <- finalize_workflow(wf, best_params)
    final_fit <- fit(final_wf, data = train)

    preds <- predict(final_fit, new_data = test) %>%
      dplyr::pull(.pred)

    return(preds)
  }

  # Run all four ML models
  pred_rf <- tune_predict(rf_spec, "RF")
  pred_xgb <- tune_predict(xgb_spec, "XGB")
  pred_knn <- tune_predict(knn_spec, "KNN")
  pred_svm <- tune_predict(svm_spec, "SVM")

  # Stop parallel processing
  stopCluster(cl)
  registerDoSEQ()

  # Add predictions to test data and convert daily dose to the administered dose
  test <- test %>%
    dplyr::mutate(
      RF = exp(pred_rf) * (II / 24),
      XGB = exp(pred_xgb) * (II / 24),
      KNN = exp(pred_knn) * (II / 24),
      SVM = exp(pred_svm) * (II / 24),
      DOSE_ADM = DOSE_ADM *  (II / 24),
      DOSE_TARGET =  DOSE_TARGET *  (II / 24)
    )

  return(test)
}
