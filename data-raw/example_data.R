# source("data-raw/example_data.R")
# Code to generate the data used for the function examples 
library(stars)
random <- function(cont = TRUE) {
  m <- secr::make.mask(nx = 10, ny = 10, spacing = 2)
  h <- secr::randomHabitat(m, p = 0.5, A = 0.3)
  r <- secr::raster(h)

  if (cont) {
    r <- r * rnorm(length(r), 1, .2)
    r <- r / raster::maxValue(r)
  }
  
  stars::st_as_stars(r)
}

# Drivers 
nDr <- 10
drivers <- purrr::map(1:nDr, \(i) random()) |>
           setNames(glue::glue("driver{1:nDr}"))
drivers <- do.call("c", drivers) 
# plot(merge(drivers))

# Valued components
nVC <- 20
vc <- purrr::map(1:nVC, \(i) random(FALSE)) |>
           setNames(glue::glue("vc{1:nVC}"))
vc <- do.call("c", vc) 
# plot(merge(vc))

# Sensitivity 
sensitivity <- matrix(
  ncol = nVC,
  nrow = nDr,
  data = runif(n = nVC*nDr),
  dimnames = list(names(drivers), names(vc))
)

# Metaweb 
metaweb <- matrix(
  ncol = nVC,
  nrow = nVC,
  data = sample(c(0,1), size = nVC * nVC, replace = TRUE),
  dimnames = list(names(vc), names(vc))
)
diag(metaweb) <- 0

# Export
usethis::use_data(drivers, vc, sensitivity, metaweb, overwrite = TRUE, internal = TRUE)
