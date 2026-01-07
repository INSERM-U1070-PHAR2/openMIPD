#' Compute probability of target attainment based on fToverMIC.
#'
#' @param tbl_t_over_mic Tibble with columns ft_over_mic and MIC.
#' @param target Target fraction of time above MIC.
#'
#' @return Tibble with fraction of target attainment for each MIC value in
#'   tbl_t_over_mic.
#' @export
#' @importFrom dplyr summarise
#' @examples
#' compute_pta_t_above_mic(tbl_t_over_mic = tbl_t_over_mic, target = 1)

compute_pta_t_above_mic <- function(tbl_t_over_mic,
                                   target) {
  tbl_t_over_mic |>
    dplyr::summarise(
      .by = MIC,
      mean_ft_over_mic = mean(ft_over_mic),
      sd_ft_over_mic = sd(ft_over_mic),
      median_ft_over_mic = median(ft_over_mic),
      "PTA_{target}" := sum(ft_over_mic >= target) / dplyr::n()
    )
}
