#' Execute a complete assessment workflow
#'
#' This function is used execute an entire assessment workflow including the evaluation of stressors and valued component cumulative footprint, the cumulative exposure of valued components to stressors, and the type cumulative effects assessment selected (i.e. individual or network). It uses an assessment workflow yaml configuration file to execute the complete the workflow.
#'
#' @param config path to a yaml assessments workflow configuration file prepared by the user. Use `new_assessment()` to generate a new configuration file template. Alternatively, this can be a list organized as a yaml document.
#'
#' @return This function returns data and figures from the assessment with all the relevant elements specified in the yaml configuration file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' assessment(config = "./assessment.yml")
#' }

assessment <- function(config) {
  # Load yaml configuration file 
  dat <- yaml::read_yaml(config)

  # Stressors 
  stressors <- dat$assessment$stressors
  nst <- length(stressors)
  stid <- character(nst)
  for(i in 1:nst) stid[i] <- stressors[[i]]$id
  dat <- pipedat::importdat(stid)
  
  vc <- dat$assessment$valued_components
  vc[[1]][[1]]$id
  
  for(i in 1:length(stressors)) 
  stressors[[1]]
  
  
  
  # Parameters
  crs <- dat$data_workflow$parameters$crs
  bbox <- unlist(dat$data_workflow$parameters$bbox)
  timespan <- dat$data_workflow$parameters$timespan
  params <- dat$data_workflow$params

}
