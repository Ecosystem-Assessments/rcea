#' Cumulative footprint
#'
#' Function to evaluate cumulative footprint, i.e. the sum of individual elements or their intensity at a given location, of stressors or valued components
#'
#' @param dat data.frame with stressors intensity or valued components, with two columns: `uid` for the unique identifier of the grid cells; `value` for the value of the stressors or valued components that should be cumulated
#'
#' @export
#'
#' @example 
#' dat <- data.frame(uid = rep(1:5,2), value = runif(10, 0, 100))
#' cumulfoot(dat)

cumulfoot <- function(dat) {
  dplyr::group_by(dat, uid) |>
  dplyr::summarise(cumulative_footprint = sum(value, na.rm = TRUE)) |>
  data.frame()
}
