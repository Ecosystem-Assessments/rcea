#' `rcea` helper functions
#'
#' List of functions that support the main functions of the `rcea` package
#'
#' @param dat list of cea matrices or cea array
#' @eval arguments(c("drivers","vc"))
#'
#' @describeIn helpers create array from list of cea matrices
#' @export
make_array <- function(dat) {
  unlist(dat) |>
  array(
    dim = c(
      nrow(dat[[1]]), 
      ncol(dat[[1]]), 
      length(dat)
    ),
    dimnames = list(
      c(), 
      names(dat[[1]]), 
      names(dat)
    )
  )  
}

#' ========================================================================================
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' ----------------------------------------------------------------------------------------
#' @describeIn helpers create stars object from list of cea matrices or cea array
#' @export
make_stars <- function(dat, drivers, vc) {
  xy <- sf::st_coordinates(drivers)
  drNames <- names(drivers)
  vc_index <- data.frame(
    vc = names(vc),
    vc_id = 1:length(vc)
  )
    
  # For data exported as list format (exposure and cea assessments)
  if ("list" %in% class(dat)) {
    for(i in 1:length(dat)) {
      dat[[i]] <- cbind(xy, dat[[i]]) |>
                  dplyr::mutate(drivers = drNames[i])
    }
    dat <- dplyr::bind_rows(dat) |>
           stars::st_as_stars(dims = c("x","y","drivers"))          
  }
  
  # For data exported in long format, from the ncea assessment
  if ("data.frame" %in% class(dat)) {
    # Get coordinates with repeated driver names
    xy$id_cell <- 1:nrow(xy)
    xyd <- cbind(
      dplyr::slice(xy, rep(1:dplyr::n(), each = length(drNames))),
      drivers = rep(drNames, nrow(xy))
    )
    
    # Pivot data to get drivers as lines and add vc names
    dat <- tidyr::pivot_longer(
      dat,
      cols = -c(id_cell, vc_id), 
      names_to = "drivers", 
      values_to = "effect"
    ) |>
    dplyr::left_join(vc_index, by = "vc_id") |>
    dplyr::select(-vc_id) |>
    tidyr::pivot_wider(names_from = vc, values_from = effect) 
    
    # Join with xy data and transform into stars object
    dat <- dplyr::left_join(xyd, dat, by = c("id_cell", "drivers")) |>
           dplyr::select(-id_cell) |>
           dplyr::select(x, y, drivers, dplyr::all_of(vc_index$vc)) |>
           stars::st_as_stars(dims = c("x","y","drivers"))          
  }
  
  # Return 
  dat
}