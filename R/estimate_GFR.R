#' Estimate glomerular filtration rate (GFR) using the 2021 CKD-EPI equation.
#'
#' @param AGE Age in years.
#' @param CREAT Serum creatinine concentration in mg/dL.
#' @param SEX "M" or "F".
#'
#' @returns Estimated GFR values in mL/min/1.73 m².
#' @export
#' @examples
#' data$eGFR <- estimate_GFR(data$AGE, data$CREAT, data$SEX)

estimate_GFR <- function(AGE, CREAT, SEX) {
  if (SEX %in% c("M", "F")) {

    # Set A and B based sex and the value of serum creatinine
    if (SEX == "F") {
      if (CREAT <= 0.7) {
        A <- 0.7
        B <- -0.241
      } else {
        A <- 0.7
        B <- -1.2
      }
    } else if (SEX == "M") {
      if (CREAT <= 0.9) {
        A <- 0.9
        B <- -0.302
      } else {
        A <- 0.9
        B <- -1.2
      }
    }

    # Mutliply by 1.012 for females
    sex_factor <- ifelse(SEX == "M", 1, 1.012)

    # Calculate eGFR
    eGFR <- 142 * ((CREAT / A)^B) * (0.9938^AGE) * sex_factor
  } else {
    eGFR <- NA
  }

  return(eGFR)
}
