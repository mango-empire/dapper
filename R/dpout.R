#' Creates `dpout` object.
#'
#' @param theta parameter being estimated.
#' @param accept_mat acceptance matrix.
#' @param varnames variable names.
#'
#' @return a simulation object
#' @noRd
new_dpout <- function(theta, accept_mat, varnames = NULL) {
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
  e1 <- posterior::as_draws_matrix(dp_obj)
  e2 <- accept_mat
  structure(list(chain = e1, accept_prob = accept_mat), class = c("dpout"))
}


