# ------------------------------------------------------------------------------
# List of arguments used throughout the package
# NOTE: Documented here in order to avoid unnecessary repetition in documentation and allow for easy modifications in the future if necessary.
arguments <- function(arg) {
  dat <- list()
  dat$drivers <- "@param drivers list of stars or terra objects"
  dat$vc <- "@param vc list of stars or terra objects"
  dat$sensitivity <- "@param sensitivity matrix of drivers x vc"
  dat$metaweb <- "@param metaweb matrix of vc x vc"
  dat$trophic_sensitivity <- "@param trophic_sensitivity data.frame of trophic sensitivities, default from Beauchesne"
  dat$weight <- "@param weight weight for the direct vs indirect modules when calculating cea scores, should be between 0 and 1"
  
  dat[names(dat) %in% arg] |>
  unlist()
}
