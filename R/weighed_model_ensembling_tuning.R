#' Automatic hyperparameter (RMSE and MPE penalty) tuning for weighed model ensembling. Selection criterion is target attainment in the training data.
#'
#' @param train_data Training dataset with ID to indicate a subject (set of covariates), covariates, observed target variable (for example CMAX), predicted target variable (for example CMAX_PRED), TIME if the target_variable is CONCENTRATION.
#' @param target_variable A secondary PK parameter (CMAX - maximal concentration, AUC - exposure or CMIN - trough concentration) or CONCENTRATION.
#' @param continuous_cov List of continuous covariates.
#' @param categorical_cov List of categorical covariates.
#' @param penalties_grid Tibble with two columns pen_RMSE and pen_MPE containing penalties values to be tested.
#' @param conc_inf lower limit of target concentration interval.
#' @param conc_sup upper limit of target concentration interval.
#' @param conc_target target concentration.
#' @param freq_column name of the interdose interval column.
#' @param cv_stratification_col name of the column indicating the stratification levels for stratified cross-validation. Defaults to NULL.
#' @param n_cv number of cross validation folds. Defaults to 10.
#' @param target_log_dose should the tuning be based on log dose or dose ? Default to TRUE.

#' @returns Tuned rRMSE and MPE penalties.
#' @export
#' @import purrr
#' @import dplyr
#' @examples
#' best <- weighed_model_ensembling_tuning(train_data = AMOX_CMIN_TRAIN, test_data = AMOX_CMIN_TEST, target_variable = "CMIN", continuous_cov = c("WT", "CREAT", "AGE"),categorical_cov = c("OBESE", "ICU", "BURN", "SEX"), conc_inf = 40, conc_sup = 80, conc_target = 60, max_dose_g = 20, freq_column = "FREQ")

weighed_model_ensembling_tuning <- function(
    train_data,
    target_variable,
    continuous_cov, categorical_cov,
    penalties_grid = NULL,
    conc_inf = 40, conc_sup = 80, conc_target = 60,
    freq_column = "FREQ",
    cv_stratification_col = NULL,
    n_cv = 10,
    target_log_dose = TRUE) {

  if(is.null(penalties_grid)){
  # Penalty combinations to test
  penalties_grid <- tibble::tribble(
    ~pen_RMSE, ~pen_MPE,
    -4, -20,
    -25, -5,
    -25, -3,
    -25, -1,
    -4, -2,
    -8, -1,
    -2, -4,
    -2, -6,
    -2, -8,
    -1, -8
  )}

  folds <- rsample::vfold_cv(train_data, v = n_cv, strata = cv_stratification_col)
  pen_RMSE <- penalties_grid$pen_RMSE
  pen_MPE <- penalties_grid$pen_MPE
  cv_train_data_list <- map(folds$splits,rsample::analysis)
  cv_test_data_list <- map(folds$splits,rsample::assessment)

  test_results <- pmap( # pmap to make it faster than a loop
    list(rep(pen_RMSE,each=n_cv),
         rep(pen_MPE,each=n_cv),
         rep(cv_train_data_list,n_cv),
         rep(cv_test_data_list,n_cv)),
    function(pen_RMSE, pen_MPE,cv_train_data,cv_test_data) {

      # WME train
      train_out <- weighed_model_ensembling_train(
        data = cv_train_data,
        target_variable = target_variable,
        pen_RMSE = pen_RMSE,
        pen_MPE = pen_MPE,
        continuous_cov = continuous_cov,
        categorical_cov = categorical_cov
      )

      # WME test applied to training data
      test_out <- weighed_model_ensembling_test(
        test_data = cv_test_data,
        train_results = train_out
      )

      # Calculate target attainment for training data
      test_eval <- test_out$test_results |>
       dplyr::filter(REFERENCE == 1) |>
        mutate(
          DOSE_PRED = (60 / WEIGHTED_PREDICTION) * DOSE_ADM,
          DOSE_TRUE = (60 / CMIN_IND) * DOSE_ADM,
          DOSE_inf  = (conc_inf / CMIN_IND) * DOSE_ADM,
          DOSE_sup  = (conc_sup / CMIN_IND) * DOSE_ADM,
          log_DOSE_TRUE = log(DOSE_TRUE),
          log_DOSE_PRED = log(DOSE_PRED),
          Prediction_correctness = ifelse(
            (DOSE_PRED >= DOSE_inf & DOSE_PRED <= DOSE_sup),
            "Correct", "Incorrect"
          )
        )
if(target_log_dose == TRUE){
      rmse_test <- yardstick::rmse(test_eval,
                                   truth = log_DOSE_TRUE,
                                   estimate = log_DOSE_PRED)
} else
{
  rmse_test <- yardstick::rmse(test_eval,
                               truth = DOSE_TRUE,
                               estimate = DOSE_PRED)
}

      attainment <- mean(test_eval$Prediction_correctness == "Correct",
                         na.rm = TRUE)

      list(
        train_out = train_out,
        test_out = test_eval,
        target_attainment = attainment,
        pen_MPE = pen_MPE,
        pen_RMSE = pen_RMSE,
        rmse_test = rmse_test
      )
    }
  )

  summary_test_results <- tibble(correct = map(test_results,~.x$target_attainment) |> unlist(),
                                 rmse = map(test_results,~.x$rmse_test$.estimate) |> unlist(),
                         pen_MPE = map(test_results,~.x$pen_MPE) |> unlist(),
                         pen_RMSE = map(test_results,~.x$pen_RMSE) |> unlist(),
                         fold = rep(1:n_cv,n_cv)) |>
    group_by(pen_MPE,pen_RMSE) |>
    summarise(mean_precision = mean(correct),
              sd_precision = sd(correct),
              mean_rmse = mean(rmse),
              sd_rmse = sd(rmse),
              median_rmse = median(rmse)) |>
    arrange(median_rmse) |>
    ungroup()

  # Select the best penalty combination
  best_pens <- summary_test_results |>
   dplyr::filter(median_rmse == min(median_rmse))

  if(target_log_dose == TRUE){
    message("Selected penalties: pen_RMSE = ", best_pens$pen_RMSE ,
            ", pen_MPE = ", best_pens$pen_MPE,
            " (Mean target attainment = ", round(best_pens$mean_precision*100, 1), "%)",
            " (Median RMSE (prediction of log(DOSE)) = ", best_pens$median_rmse, ")")
  } else
  {
    message("Selected penalties: pen_RMSE = ", best_pens$pen_RMSE ,
            ", pen_MPE = ", best_pens$pen_MPE,
            " (Mean target attainment = ", round(best_pens$mean_precision*100, 1), "%)",
            " (Median RMSE (prediction of DOSE) = ", best_pens$median_rmse, ")")
  }



  return(
    list(cv_results = summary_test_results,
         results_best_pens = best_pens,
         best_pen_MPE = best_pens$pen_MPE,
         best_pen_RMSE = best_pens$pen_RMSE)
    )
}

