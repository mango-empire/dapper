.onLoad <- function(libname, pkgname) {
    ddnorm_constant <<- memoise::memoise(ddnorm_constant)
}
