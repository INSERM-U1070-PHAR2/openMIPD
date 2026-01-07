#' XGBoost (XGB) for MIPD
#'
#' @param final_fit The trained XGBoost model from the xgb_train function.
#' @param test A preprocessed (ml_data_preprocess function) training dataset with one line per patient.
#' @param seed Seed for reproducibility.
#' @returns A test data with containing the XGBoost predictions.
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
#' @import stacks
#' @examples
#' results <- xgb_test(final_fit = final_fit, test = AMOX_CMIN_TEST)

xgb_test <- function(final_fit,
                test,
                seed = 1991) {

  set.seed(seed)

preds <- predict(final_fit, new_data = test) %>%
  dplyr::pull(.pred)

# Collect predictions into test set
test_final <- test %>%
  dplyr::mutate(
    XGBoost = exp(preds) * (II / 24),
    DOSE_ADM = DOSE_ADM *  (II / 24),
    DOSE_TARGET = DOSE_TARGET *  (II / 24)
  )

return(test_results = test_final)

}
