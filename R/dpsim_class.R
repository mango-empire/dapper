#' Creates simulation object
#'
#' @param theta
#'
#' @return a simulation object
#' @export
#'
#' @examples
new_dpsim <- function(theta, varnames = NULL) {
  dp_obj <- do.call(rbind, theta)
  nr <- length(theta) * nrow(theta[[1]])
  nc <- ncol(theta[[1]])
  vn <- paste0("theta", 1:nc)
  if(!is.null(varnames)) vn <- varnames
  dl <- list(draw = as.character(1:nr),
             variable = vn)
  attr(dp_obj, 'dimnames') <- dl
  attr(dp_obj, 'nchains') <- length(theta)
  #structure(dp_obj, class = c("dpsim", "matrix"))
  posterior::as_draws_matrix(dp_obj)
}


