#' Data preprocessing for machine learning
#'
#' @param target_variable A secondary PK parameter (CMAX - maximal concentration, AUC - exposure or CMIN - trough concentration).
#' @param data A dataset with one line per patient. The dataset should include the target variable (ending with _IND for example CMIN_IND), and the administered dose dose (DOSE_ADM) and a column called II for interdose interval (set to 24 for continuous administration).
#' @param target_concentration Target to reach (Cmax, Cmin or AUC).
#' @returns A test data the log transformed target dose, and dosing scheme based stratification.
#' @export
#' @import dplyr
#' @examples
#' results <- ml_data_preprocess(data = AMOX_CMIN_TRAIN, target_variable = "CMIN", target_concentration = 60)

ml_data_preprocess <- function(data,
                             target_variable,
                             target_concentration) {


  target_variable <- paste0(target_variable, "_IND") # Individual concentrations end with _IND

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

  data_preprocessed <- preprocess_data(data)
  return(data_preprocessed)
}
