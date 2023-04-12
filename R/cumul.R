#' Cumulative footprint
#'
#' Function to evaluate cumulative footprint, i.e. the sum of the distribution and intensity/probability of stressors or valued components at a given location
#'
#' @param dat list of stars or terra objects to measure cumulative footprints
#'
#' @export
cumul <- function(dat) {
  dplyr::group_by(dat, uid) |>
  dplyr::summarise(cumulative_footprint = sum(value, na.rm = TRUE)) |>
  data.frame()
}
