# glm4

## Overview

`glm4` is an R package for fitting generalised linear models (GLMs) on large, high-dimensional datasets where standard approaches run out of memory.

This package adds functionality to `MatrixModels::glm4()`, which does the core work of fiting GLMs using iteratively reweighted least squares (IRLS) with support for sparse matrix arithmetic via the `Matrix` package. When a model includes many factor levels (for example, thousands of fixed effects) the usual dense format model matrix can require tens or hundreds of gigabytes of RAM. Sparse storage keeps only the non-zero elements, making fitting and working with these models more feasible.

The raw output of `MatrixModels::glm4()` is a bare-bones S4 object with no print, summary, or inference methods. `glm4` fills that gap by wrapping the fit in a familiar S3 object with the same structure and methods as a standard `stats::glm` fit, including:

- `print()` and `summary()` — coefficients, standard errors, and test statistics
- `logLik()`, `AIC()`
- `confint()` — Wald confidence intervals
- `vcov()` — variance-covariance matrix (returned as a sparse `Matrix`)
- `predict()`, `residuals()`
- `anova()` — sequential deviance analysis and multi-model comparison, matching `stats::anova.glm` output exactly

The interface is intentionally close to `stats::glm()`. In most cases, switching from `glm()` to `glm4()` requires only adding `sparse = TRUE`.

## Installation

`glm4` is not yet on CRAN. Install the development version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("awhug/glm4")
```

## Contributing

Contributions are more than welcome. If you find a bug or have a feature request, please [open an issue](https://github.com/awhug/glm4/issues). If you would like to contribute code, feel free to submit a pull request.
