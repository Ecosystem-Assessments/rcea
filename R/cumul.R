#' Cumulative footprint
#'
#' Function to evaluate cumulative footprint, i.e. the sum of the distribution and intensity/probability of stressors or valued components at a given location
#'
#' @param dat list of stars or terra objects to measure cumulative footprints
#'
#' @export
cumul <- function(dat) {
  library(stars)
  do.call("c", dat) |>
  stars::st_redimension() |>
  stars::st_apply(c(1,2), sum, na.rm = TRUE)      
}
