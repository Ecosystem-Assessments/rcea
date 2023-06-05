#' rcea: package to perform cumulative effects assessments
#'
#' @docType package
#' @name rcea
#'
#' @importFrom glue glue glue_sql
#' @importFrom grDevices dev.off png colorRampPalette
#' @importFrom graphics par box layout lines mtext
#' @importFrom graphics polygon text
#' @importFrom rlang sym
#' @importFrom yaml yaml.load_file write_yaml read_yaml
#' @importFrom cli symbol
#' @importFrom crayon blue
NULL


# ------------------------------------------------------------------------------
# Gracieuseté de Kevin Cazelles: https://github.com/KevCaz
# my simple(r) version of use template
use_template <- function(template, save_as = stdout(), pkg = "rcea", ...) {
  template <- readLines(
    fs::path_package(package = pkg, template)
  )
  # NB by default whisker forward the parent envi and I used this
  writeLines(whisker::whisker.render(template, ...), save_as)
}

#' Check if folder exists and create if not
#'
#' @param path path of folder to use as output, create if it does not already exist
#'
#' @export
chk_create <- function(path) {
  if (!file.exists(path)) dir.create(path, recursive = TRUE)
}