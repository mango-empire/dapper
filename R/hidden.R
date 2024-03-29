
parse_formula <- function(formula, data=NULL, intercept=TRUE) {

  # extract Y, X, and variable names for model formula and frame
  mf <- match.call(expand.dots = FALSE)
  mf$intercept <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  mt <- attr(mf, "terms")
  if (!intercept){
    attributes(mt)$intercept <- 0
  }

  # null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf)
  X <- as.matrix(X)         # X matrix
  xvars <- dimnames(X)[[2]] # X variable names
  xobs  <- dimnames(X)[[1]] # X observation names
  Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
  list(Y, X, xvars, xobs)
}
