#' Calculate body surface area (BSA) based on the Du Bois & Du Bois equation.
#'
#' @param WT Body weight in kg.
#' @param HT Height in cm.
#'
#' @returns Calculated BSA in m².
#' @export
#' @examples
#' data$BSA <- calculate_BSA(data$WT, data$HT)

calculate_BSA <- function(WT, HT) {

    # Calculate eGFR
    BSA <- ((WT^0.425) * (HT^0.725) * 0.007184)


  return(BSA)
}
