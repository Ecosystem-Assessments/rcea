#' Creates new files from templates
#'
#' This function is used to create new files from templates available in the `rcea` package.
#'
#' @param name name of the data workflow; defaults to `assessment`
#'
#' @return This function creates new files from the template available in `inst/templates/`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # New assessment workflow
#' new_assessment()
#' }
new_assessment <- function(name = "assessment") {
  use_template(
    template = "templates/assessment.yml",
    save_as = glue::glue("./{name}.yml")
  )
}
