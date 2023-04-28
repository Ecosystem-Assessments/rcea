#' Cumulative effects assessments
#'
#' Assessment of cumulative effects using the Halpern et al. 2008 method. 
#'
#' @eval arguments(c("drivers","vc","sensitivity"))
#' @param exportAs string, the type of object that should be created, either a "list" or a "stars" object
#'
#' @export
#'
#' @examples
#' # Data
#' drivers <- rcea:::drivers 
#' vc <- rcea:::vc
#' sensitivity <- rcea:::sensitivity
#'
#' # Species-scale effects
#' (halpern <- cea(drivers, vc, sensitivity, "stars"))
#' plot(halpern)
#' halpern <- merge(halpern, name = "vc") |>
#'           split("drivers")
#' plot(halpern)
cea <- function(drivers, vc, sensitivity, exportAs = "list") {
  # Exposure
  dat <- exposure(drivers, vc)

  # Sensitivity
  nmDr <- names(dat)
  nmVC <- names(dat[[1]])
  sensitivity <- sensitivity[nmVC, nmDr] 
  
  # Effect of drivers on valued components (D * VC * u)
  for(i in 1:length(dat)) {
    dat[[i]] <- sweep(dat[[i]], MARGIN = 2, sensitivity[,i], `*`)
  }

  # Return
  if (exportAs == "list") {
    dat
  } else if (exportAs == "stars") {
    xy <- sf::st_coordinates(vc)
    drNames <- names(drivers)
    make_stars(dat, drivers, vc)
  }

}