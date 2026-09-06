# dapper
<!-- badges: start -->
[![codecov](https://codecov.io/gh/mango-empire/dapper/graph/badge.svg?token=GP7KF5RR3O)](https://app.codecov.io/gh/mango-empire/dapper)
<!-- badges: end -->

A data augmentation framework for privacy aware Bayesian inference.

## Installation
Install the development version, including its vignette, with:

``` r
devtools::install_github("mango-empire/dapper", build_vignettes = TRUE)
```

## Example workflow

The randomized-response vignette walks through model construction, posterior
sampling, diagnostics, and a comparison with an analysis that ignores privacy
noise. After installation, open it with:

``` r
vignette("randomized-response", package = "dapper")
```
