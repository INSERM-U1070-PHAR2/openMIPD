#' Plot to show average model weights for categories/quantiles of a covariate in the test data
#'
#' @param selected_regiments List of regimen_ids (dosing regimen identifiers)to use in the simulation.
#' @param selected_IDs_dose List of which IDs (subjects) to include in the simulation.
#' @param dosing_table A dataset with the characteristics of the dosing regimens containing the following columns: regimen_id (dosing regimen identifier), time, tinf (duration of infusion), amt, ii (interdose interval), evid, cmt (compartment).
#' @param COV_table Is the covariate categorical? TRUE or FALSE, default is FALSE
#' @param observation_grid A list of time values for the observations.
#' @param model_list List of models to use for the simulation. The models' .cpp file should be provided by the user.
#' @returns A dataset with PRED, IPRED & Y
#' @export
#' @import dplyr
#' @import tidyr
#' @import mrgsolve
#' @import purrr
#' @examples
#' simulated_set <- simulate_mrgsolve(selected_regimens = c(10000mg, 20000mg), selected_IDs_dose = c(1:2500), dosing_table = dosing_table, COV_table = Cov_table, observation_grid = seq(0, 24, by = 2), model_list = (model1, model2))

simulate_mrgsolve <- function(selected_regimens, selected_IDs_dose, dosing_table, COV_table,

                              observation_grid, model_list) {

  results <- list()

  for (regimen in selected_regimens) {
    dosing_table_filt <- dosing_table |> dplyr::filter(regimen_id == regimen)

    dosing_grid1 <- dosing_table_filt |>
      group_by(tinf, amt, ii, evid, cmt) |>
      summarise(time = seq(0, 24 - ii, by = ii), .groups = "drop") |>
      mutate(rate = amt / tinf) |>
      select(time, amt, rate, ii, evid, cmt)

    final_dosing_grid <- crossing(ID = selected_IDs_dose, dosing_grid1) |>
      full_join(COV_table, by = "ID")

    event_data <- bind_rows(final_dosing_grid, observation_grid) |>
      arrange(ID, time, desc(evid))

    simulate_model <- function(model, model_name, is_zero_re) {
      model |>
        data_set(event_data) |>
        mrgsim() |>
        as.data.frame() |>
        dplyr::filter(time == 24) |>
        dplyr::select(ID, TIME = time, IPRED) |>
        dplyr::rename(PRED = IPRED) |>
        dplyr::mutate(MODEL = model_name, zero_re = is_zero_re, regimen = regimen)
    }

    # Run simulations for each model in the provided list
    data <- purrr::map_dfr(names(model_list), function(model_name) {
      models <- model_list[[model_name]]
      bind_rows(
        simulate_model(models$model, model_name, FALSE),
        simulate_model(models$zero_re, model_name, TRUE)
      )
    })

    results[[regimen]] <- data
  }

  # Combine results and merge covariates
  final_results <- bind_rows(results) |>
    left_join(COV_table, by = "ID")

  return(final_results)
}
