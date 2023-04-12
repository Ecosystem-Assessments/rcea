#' Cumulative footprint
#'
#' Function to evaluate cumulative footprint, i.e. the sum of the distribution and intensity/probability of stressors or valued components at a given location
#'
#' @param dat list of stars or terra objects to measure cumulative footprints
#'
#' @export
cumul <- function(dat) {
  library(stars)
  
  # Combine if in list
  if (class(dat) == "list") dat <- do.call("c", dat)

  # Redimension and sum 
  stars::st_redimension(dat) |>
  stars::st_apply(c(1,2), sum, na.rm = TRUE) |>
  setNames("Footprint")     
}
