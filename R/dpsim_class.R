#' Creates simulation object
#'
#' @param theta
#'
#' @return a simulation object
#' @export
#'
#' @examples
new_dpsim <- function(theta) {
  plist <- list(theta = theta)
  structure(plist, class = "dpsim")
}
