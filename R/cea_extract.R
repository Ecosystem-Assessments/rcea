#' Extract attributes and aggregate cumulative effects assessment results
#'
#' Function to extract specific drivers and vc of interest from the drivers, valued component, exposure and cumulative effects assessments results, and to aggregate data over drivers, valued components, or drivers and valued components.
#'
#' @param dat stars object, either the drivers, valued components, exposure, cumulative effects assessment or network-sacale cumulative effects assessment results
#' @param dr_sel string, name of drivers to extract 
#' @param vc_sel string, name of valued components to extract
#' @param cumul_fun function to apply on the stars object, one of "drivers" for the cumulative effects of each drivers on all value components, "vc" for the cumulative effects of all drivers on each valued component, "full" for the cumulative effects of all drivers on all valued components, "footprint" for the cumulative footprint of drivers or valued components, and "none" to keep data as is, i.e. to extract specific drivers and valued components from the assessment results.
#'
#' @examples
#' # Data
#' drivers <- rcea:::drivers 
#' vc <- rcea:::vc
#' sensitivity <- rcea:::sensitivity
#' 
#' # Cumulative effects assessment
#' dat <- cea(drivers, vc, sensitivity, "stars")
#' 
#' # Extract attributes
#' dr_sel <- c("driver1","driver5")
#' vc_sel <- c("vc4","vc7","vc10","vc12")
#' cea_extract(dat, dr_sel = dr_sel, vc_sel = vc_sel) 
#'
#' # Cumulative footprint of selected drivers and valued components
#' cea_extract(drivers, dr_sel = dr_sel, cumul_fun = "footprint") 
#' cea_extract(vc, vc_sel = vc_sel, cumul_fun = "footprint") 
#'
#' # Cumulative effects of all drivers on all vc
#' cea_extract(dat, cumul_fun = "drivers")
#'
#' # Cumulative effects of all drivers on each vc
#' cea_extract(dat, cumul_fun = "vc") 
#'
#' # Full cumulative effects
#' cea_extract(dat, cumul_fun = "full") 
#'
#' @export
cea_extract <- function(dat, dr_sel = NULL, vc_sel = NULL, cumul_fun = "none") {
  # Select drivers and vc to extract
  dat <- select_attr(dat, dr_sel, vc_sel)
  
  # Apply relevant aggregations, if applicable
  dat <- aggr(dat, cumul_fun)
  
  # Return
  dat
}

select_attr <- function(dat, dr_sel = NULL, vc_sel = NULL) {
  # Select drivers, if applicable
  if (!is.null(dr_sel)) {
    if ("drivers" %in% names(stars::st_dimensions(dat))) {
      dat <- dat[,,,dr_sel]      
    } else {
      dat <- dat[dr_sel]
    }
  }
  
  # Select valued components, if applicable
  if (!is.null(vc_sel)) dat <- dat[vc_sel]
  
  # Return
  dat
}

cumul <- function(dat) {
  stars::st_apply(dat, c("x","y"), sum, na.rm = TRUE)
}

cumul_vc <- function(dat) {
  cumul(dat)
}

cumul_drivers <- function(dat) {
  merge(dat, name = "vc") |>
  split("drivers") |>
  cumul()
}

cumul_full <- function(dat) {
  dat |>
  cumul() |>
  merge() |>
  cumul() |>
  setNames("cumulative_effects")
}

cumul_footprint <- function(dat) {
  merge(dat) |>
  cumul() |>
  setNames("cumulative_footprint")
}

aggr <- function(dat, cumul_fun) {
  switch(
   cumul_fun, 
   "drivers" = cumul_drivers(dat),
   "vc" = cumul_vc(dat),
   "full" = cumul_full(dat),
   "footprint" = cumul_footprint(dat),
   "none" = dat
  )
}
