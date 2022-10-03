#' Data normalization
#'
#' Function to scale data between 0 and 1 using the 99th percentile
#'
#' @param x numeric vector of data to normalize 
#'
#' @rdname data_transformation
#'
#' @export
#'
#' @examples
#' dat <- runif(10, 0, 100)
#' x <- quantNorm(dat)
#' cbind(dat,x)

quantNorm <- function(x) {
  id <- x != 0
  x <- x / quantile(x[id], probs = .99, na.rm = T)
  x[x > 1] <- 1
  x[x < 0] <- 0
  x
}

#' Data transformation
#'
#' Wrapper function to log transform data. This function could/should become more complex eventually
#'
#' @param x numeric vector of data to transform 
#'
#' @rdname data_transformation
#'
#' @export
#'
#' @examples
#' dat <- runif(10, 0, 100)
#' x <- logTrans(dat)
#' cbind(dat,x)

logTrans <- function(x) {
  log(x + 1)
}

#' Data normalization and transformation
#'
#' Function to normalize and transform data for further analyses
#'
#' @param x numeric vector of data to transform 
#'
#' @rdname data_transformation
#'
#' @export
#'
#' @examples
#' dat <- runif(10, 0, 100)
#' x <- dataTrans(dat)
#' cbind(dat,x)

dataTrans <- function(x) {
  logTrans(x) |>
  quantNorm()
}



